# DACSIF ALPlotter Web

ALPlotter Web is an interactive chemoinformatics web application developed using Python and Streamlit. Its main objective is to facilitate the analysis of structure-activity relationships (SAR), the exploration of activity landscapes, and the identification of *activity cliffs* through the exhaustive calculation of paired comparisons and the SALI index.

## ✨ Key Features

* **Multiple Fingerprint Calculation:** Native support via RDKit for Morgan (Circular) fingerprints, 2D Pharmacophore, MACCS Keys, Topological (RDKit), and Atom Pairs.
* **Dimensionality Reduction (t-SNE):** Visualization of the global chemical space projected in 2D, with interactivity to inspect identifiers (IDs) without overloading the browser’s memory.
* **SAS (Structure-Activity Similarity) Map:** An interactive graph that correlates structural similarity (Tanimoto) with differences in biological activity.
* **Visual Analysis by SALI Index:** A paired inspection module that automatically ranks the most relevant *Activity Cliffs* and renders comparative 2D structures side by side.
* **Cross-Highlighting:** When you click on a pair of molecules on the SAS map, the corresponding points are automatically highlighted with special markers on the t-SNE map.
* **Automatic Data Cleaning:** Built-in SMILES code deduplication and filtering of invalid structures.
* **Data Export:** Direct download of paired results, Tanimoto indices, activity deltas, and SALI metrics in CSV format.

## 🚀 Online Deployment

This tool is optimized to run in headless cloud environments such as **Streamlit Community Cloud**.

It can be easily integrated into platforms such as Google Sites via an IFrame, ensuring that end users do not need to install local dependencies.

## 💻 Local Installation and Use

If you want to run the application on your own machine or in a development environment (such as Jupyter Notebook or VS Code), follow these steps:

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/user/alplotter-web.git](https://github.com/user/alplotter-web.git)
   cd alplotter-web

   Install the necessary dependencies:
Make sure you have an active Python virtual environment (3.8 or higher).

2. **Install the dependencies**
   ```bash
    pip install -r requirements.txt

Note for Linux: You may need to install the system’s graphics dependencies (libxrender1, libxext6, libsm6) for the RDKit drawing engine to work properly.

3. **Run the application**
   ```bash
     streamlit run app.py

Usage:
Open the provided local link (usually http://localhost:8501), upload your .csv file containing the compounds, select the columns corresponding to SMILES and the activity metric, choose your fingerprint, and run the analysis.

🛠️ Technologies Used

Python 3
Streamlit (Web user interface)
RDKit (Chemoinformatics engine and 2D rendering)
Plotly (Interactive data visualization)
Scikit-learn (t-SNE algorithm)
Pandas & NumPy (Tabular data and matrix manipulation)

Code created with the assistance of Google Gemini Pro 3.1, and is based on the original R script from DIFACQUIM-UNAM. Thanks to Dr. José Luis Medina-Franco for his recommendations.
