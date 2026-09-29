import streamlit as st
import pandas as pd
import numpy as np
import base64
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Draw
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from sklearn.manifold import TSNE
import plotly.express as px
import plotly.graph_objects as go
import io

# Inicialización de estados para interactividad cruzada
if 'analyzed' not in st.session_state:
    st.session_state.analyzed = False
if 'highlight_ids' not in st.session_state:
    st.session_state.highlight_ids = []

def reset_app():
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    st.rerun()

def mol_to_base64(mol, size=(250, 250)):
    if mol is None: return ""
    img = Draw.MolToImage(mol, size=size, wedgeBonds=True)
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode()
    return f"data:image/png;base64,{img_str}"

def get_fingerprint(mol, fp_type):
    fps = {
        "Morgan (Circular)": lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048),
        "2D Pharmacophore": lambda m: Generate.Gen2DFingerprint(m, Gobbi_Pharm2D.factory),
        "MACCS Keys": lambda m: AllChem.GetMACCSKeysFingerprint(m),
        "RDKit (Topological)": lambda m: Chem.RDKFingerprint(m),
        "Atom Pairs": lambda m: AllChem.GetHashedAtomPairFingerprintAsBitVect(m, nBits=2048)
    }
    return fps[fp_type](mol)

st.set_page_config(page_title="ALPlotter web", layout="wide")
st.title("DACSIF ALPlotter Web")

# --- Barra Lateral ---
if st.sidebar.button("🔄 Reiniciar Aplicación"):
    reset_app()

st.sidebar.markdown("---")
uploaded_file = st.sidebar.file_uploader("1. Cargar dataset (CSV)", type=["csv"])
fp_option = st.sidebar.selectbox("2. Tipo de Fingerprint:", ["Morgan (Circular)", "2D Pharmacophore", "MACCS Keys", "RDKit (Topological)", "Atom Pairs"])

if uploaded_file:
    df_raw = pd.read_csv(uploaded_file)
    smiles_col = st.sidebar.selectbox("Columna SMILES", df_raw.columns)
    activity_col = st.sidebar.selectbox("Columna Actividad", df_raw.columns)
    
    df_clean = df_raw.drop_duplicates(subset=[smiles_col]).dropna(subset=[smiles_col, activity_col]).copy()
    mols_raw = [Chem.MolFromSmiles(s) for s in df_clean[smiles_col]]
    valid_mask = [m is not None for m in mols_raw]
    df_final = df_clean[valid_mask].reset_index(drop=True)
    mols = [m for m in mols_raw if m is not None]

    if st.sidebar.button("🚀 Ejecutar Análisis") or st.session_state.analyzed:
        st.session_state.analyzed = True
        
        activities = df_final[activity_col].values
        ids_list = df_final.index.tolist()

        # --- 1. Generar t-SNE ---
        fps = [get_fingerprint(m, fp_option) for m in mols]
        tsne = TSNE(n_components=2, perplexity=min(30, len(mols)-1), random_state=42)
        coords = tsne.fit_transform(np.array([list(fp) for fp in fps]))
        
        map_df = pd.DataFrame({
            'X': coords[:, 0], 'Y': coords[:, 1],
            'Actividad': activities, 'ID': ids_list
        })
        
        # --- 2. Cálculos Pareados (SALI) ---
        res_list = []
        max_a, min_a = activities.max(), activities.min()
        rango = max_a - min_a if max_a != min_a else 1
        for i in range(len(mols)):
            for j in range(i + 1, len(mols)):
                sim_t = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                delta_a = abs(activities[i] - activities[j])
                sim_act = 1 - (delta_a / rango)
                sali = delta_a / (1 - sim_t) if sim_t < 1.0 else 0
                res_list.append({
                    'ID_A': ids_list[i], 'ID_B': ids_list[j],
                    'Tanimoto': sim_t, 'Delta_Act': delta_a, 'Sim_Actividad': sim_act, 'SALI': sali,
                    'idxA': i, 'idxB': j, 'ActA': activities[i], 'ActB': activities[j]
                })
        res_df = pd.DataFrame(res_list)
        
        # --- NUEVA SECCIÓN v7.24: Análisis Visual de Pares por Índice SALI ---
        st.markdown("---")
        st.subheader("Análisis Visual de Pares por Índice SALI")
        
        # Ordenar el DataFrame de resultados por SALI de mayor a menor
        sali_df = res_df.sort_values(by="SALI", ascending=False).reset_index(drop=True)
        
        # Crear lista descriptiva para el menú desplegable
        opciones_sali = []
        for idx, row in sali_df.iterrows():
            id_a, id_b = row['ID_A'], row['ID_B']
            sali_val = row['SALI']
            delta_act = row['Delta_Act'] 
            opciones_sali.append(f"SALI: {sali_val:.2f} | Par: {id_a} vs {id_b} | Delta Act: {delta_act:.2f}")
            
        seleccion = st.selectbox("Selecciona un par para visualizar sus estructuras 2D:", opciones_sali)
        
        if seleccion:
            # Identificar la fila correspondiente a la selección
            idx_seleccion = opciones_sali.index(seleccion)
            fila_seleccionada = sali_df.iloc[idx_seleccion]
            
            # Extraer las posiciones internas usando idxA e idxB
            idx_a, idx_b = int(fila_seleccionada['idxA']), int(fila_seleccionada['idxB'])
            id_a, id_b = fila_seleccionada['ID_A'], fila_seleccionada['ID_B']
            
            # Extraer los datos originales usando df_final
            smiles_a = df_final.iloc[idx_a][smiles_col]
            smiles_b = df_final.iloc[idx_b][smiles_col]
            act_a = df_final.iloc[idx_a][activity_col]
            act_b = df_final.iloc[idx_b][activity_col]
            
            # Generar los objetos Mol de RDKit
            mol_a = Chem.MolFromSmiles(smiles_a)
            mol_b = Chem.MolFromSmiles(smiles_b)
            
            # Renderizar imágenes de alta resolución utilizando la función base64
            img_a = mol_to_base64(mol_a, size=(350, 350))
            img_b = mol_to_base64(mol_b, size=(350, 350))
            
            # Mostrar en dos columnas
            col1, col2 = st.columns(2)
            with col1:
                st.image(img_a, use_container_width=True)
                st.markdown(f"**ID:** {id_a}")
                st.markdown(f"**Actividad:** {act_a:.2f}")
                st.caption(f"SMILES: {smiles_a}")
            with col2:
                st.image(img_b, use_container_width=True)
                st.markdown(f"**ID:** {id_b}")
                st.markdown(f"**Actividad:** {act_b:.2f}")
                st.caption(f"SMILES: {smiles_b}")

        # --- 3. Renderizado Mapa SAS (Lado Izquierdo) e Inspección (Lado Derecho) ---
        st.header("1. Mapa SAS e Inspección Estructural")
        col_map, col_details = st.columns([2, 1])

        with col_map:
            fig_sas = px.scatter(res_df, x='Tanimoto', y='Delta_Act', color='SALI', 
                                color_continuous_scale='Reds',
                                custom_data=['idxA', 'idxB', 'ID_A', 'ID_B', 'ActA', 'ActB', 'SALI'])
            fig_sas.update_traces(hovertemplate="<b>Par:</b> %{customdata[2]} - %{customdata[3]}<br><b>SALI:</b> %{marker.color:.4f}")
            selected_sas = st.plotly_chart(fig_sas, use_container_width=True, on_select="rerun")

        with col_details:
            st.subheader("Detalles del Par")
            if selected_sas and "points" in selected_sas and len(selected_sas["points"]) > 0:
                p = selected_sas["points"][0]
                idxA, idxB = p["customdata"][0], p["customdata"][1]
                idA, idB = p["customdata"][2], p["customdata"][3]
                actA, actB = p["customdata"][4], p["customdata"][5]
                sali_v = p["customdata"][6]
                
                # Actualizar resaltado para t-SNE
                st.session_state.highlight_ids = [idxA, idxB]

                st.write(f"**Análisis:** ID {idA} vs ID {idB}")
                st.write(f"**SALI:** {sali_v:.4f}")
                st.markdown("---")
                st.write(f"**ID {idA}** (Act: {actA:.2f})")
                st.image(mol_to_base64(mols[idxA]))
                st.markdown("---")
                st.write(f"**ID {idB}** (Act: {actB:.2f})")
                st.image(mol_to_base64(mols[idxB]))
            else:
                st.info("Haz clic en un punto del mapa SAS para inspeccionar el par y resaltarlo en el t-SNE.")

        # --- 4. Renderizado Mapa t-SNE (Resaltado dinámico) ---
        st.header("2. Mapa t-SNE del Espacio Químico")
        fig_tsne = px.scatter(map_df, x='X', y='Y', color='Actividad', color_continuous_scale='Viridis',
                             custom_data=['ID'])
        fig_tsne.update_traces(hovertemplate="<b>ID:</b> %{customdata[0]}")
        
        # Añadir resaltado si hay selección
        if st.session_state.highlight_ids:
            h_df = map_df.iloc[st.session_state.highlight_ids]
            fig_tsne.add_trace(go.Scatter(x=h_df['X'], y=h_df['Y'], mode='markers',
                                         marker=dict(color='red', size=18, symbol='star', 
                                                    line=dict(width=2, color='white')),
                                         name='Par Seleccionado'))
        
        st.plotly_chart(fig_tsne, use_container_width=True)

        # --- 5. Exportación ---
        st.header("3. Exportar Resultados")
        csv_buffer = io.StringIO()
        res_df.drop(columns=['idxA', 'idxB']).to_csv(csv_buffer, index=False)
        st.download_button("💾 Descargar Resultados (CSV)", csv_buffer.getvalue(), "alplotter_v7_24.csv", "text/csv")
