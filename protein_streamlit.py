import streamlit as st
import py3Dmol
from stmol import showmol
from snowflake.snowpark import Session
import requests
import pandas as pd
import json
import os
import snowflake.connector
from snowflake.snowpark.functions import col

def connection() -> snowflake.connector.SnowflakeConnection:
    if os.path.isfile("/snowflake/session/token"):
        creds = {
            'host': os.getenv('SNOWFLAKE_HOST'),
            'port': os.getenv('SNOWFLAKE_PORT'),
            'protocol': "https",
            'account': os.getenv('SNOWFLAKE_ACCOUNT'),
            'authenticator': "oauth",
            'token': open('/snowflake/session/token', 'r').read(),
            'warehouse': os.getenv('SNOWFLAKE_WAREHOUSE'),
            'database': os.getenv('SNOWFLAKE_DATABASE'),
            'schema': os.getenv('SNOWFLAKE_SCHEMA'),
            'client_session_keep_alive': True
        }
    else:
        creds = {
            'account': os.getenv('SNOWFLAKE_ACCOUNT'),
            'user': os.getenv('SNOWFLAKE_USER'),
            'password': os.getenv('SNOWFLAKE_PASSWORD'),
            'warehouse': os.getenv('SNOWFLAKE_WAREHOUSE'),
            'database': os.getenv('SNOWFLAKE_DATABASE'),
            'schema': os.getenv('SNOWFLAKE_SCHEMA'),
            'client_session_keep_alive': True
        }

    connection = snowflake.connector.connect(**creds)
    return connection

def setupsession() -> Session:
    return Session.builder.configs({"connection": connection()}).create()

# Setup a Snowflake session
session = setupsession()

def find_nth(haystack, needle, n):
    parts = haystack.split(needle, n+1)
    if len(parts)<=n+1:
        return -1
    return len(haystack)-len(parts[-1])-len(needle)

# Get protein information from uniprot
def get_desc(id):
    url1='https://rest.uniprot.org/uniprotkb/' + str(id) + '.json?fields=id,organism_name'
    data1 = json.loads(requests.get(url1).text)
    uniprotid = data1['uniProtkbId']
    orgname = data1['organism']['scientificName']
    url2='https://rest.uniprot.org/uniprotkb/' + str(id) + '?fields=cc_function&format=tsv'
    data2 = requests.get(url2).text
    fx = data2[:find_nth(data2, ".", 3)+1]
    return [uniprotid, orgname, fx]

# Visualize the predicted structure saved in PDB file
def viz(id, bck):
    with open("pdb/"+ str(id) +".pdb") as ifile:
        system = "".join([x for x in ifile])
    xyzview = py3Dmol.view(query=id)
    xyzview.addModelsAsFrames(system)
    xyzview.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
    xyzview.setBackgroundColor(bck)#('0xeeeeee')
    xyzview.spin(True)
    xyzview.zoomTo()
    showmol(xyzview, height=350, width=400) 

def get_function(seq):
    df = session.sql(f"""
                    SELECT 
                        FUNCTION
                    FROM BIONEMO_DB.PUBLIC.PROTEIN_SEQUENCE_FUNCTION
                    WHERE UNIPROTID = '{seq}'""")
    protein_function = str(df.select(col('FUNCTION')).to_pandas().values)[3:-3]
    return protein_function

# Create the Streamlit
protein_list = ['P04637','P83302']
st.sidebar.title('Show Similar Proteins')
protein = st.sidebar.text_input('Input a protein to match 3 different proteins based on the ProtT5 embeddings:')
info0 = get_desc(protein)
desc = get_function(protein)

col1, col2 = st.columns(2)
with col1:
    st.write(f"Base Protein: **{protein}**")
    viz(protein, 'white')
with col2:
    st.write(f"**UniProtId:** {info0[0]}")
    st.write(f"**Organism:** {info0[1]}")
    st.text_area(label = '**Function:**', value = desc)

    # Find similar proteins of the selected protein
df = session.sql(f"""
                    SELECT 
                        UNIPROTID, 
                        VECTOR_L2_DISTANCE(EMB,(SELECT EMB FROM BIONEMO_DB.PUBLIC.PROTEINS WHERE UNIPROTID = '{protein}')::VECTOR(FLOAT,1024)) AS DISTANCE,
                        EMB
                    FROM BIONEMO_DB.PUBLIC.PROTEINS
                    ORDER BY DISTANCE ASC
                    LIMIT 4""")
similar_protein_list = df.select(col('UNIPROTID')).to_pandas().values.tolist()
similar_protein_list_all = [x for xs in similar_protein_list for x in xs]

# Display similar proteins
cols = st.columns(3)
j = 0
for i in similar_protein_list_all:
    if i not in protein_list:
        info = get_desc(i)
        with cols[j%3]:
            st.write(f"Matching Protein: **{i}**")
            st.write(f"**UniProtId:** {info[0]}")
            st.write(f"**Organism:** {info[1]}")
            st.write(info[2])
            viz(i, 'white')
    j += 1
