# Imports
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2
from io import BytesIO
from PIL import Image
import base64
import requests
from scipy.stats import spearmanr
import io

def create_clustergram(genes, terms, scores):
    """
    Creates a clustergram based on enriched genes and their respective terms.
    
    :param genes: List of enriched genes.
    :param terms: List of terms corresponding to enriched genes.
    :param scores: List of enrichment scores for genes.
    
    :return: BytesIO object containing the clustergram.
    """
    # Create a dataframe for heatmap
    df_clustergram = pd.DataFrame({'Gene': genes, 'Term': terms, 'Score': scores})
    df_clustergram = df_clustergram.pivot("Gene", "Term", "Score")
    
    # Handle NaN or infinite values
    df_clustergram = df_clustergram.fillna(0)  # fill NaN with zeros
    df_clustergram = df_clustergram.replace([np.inf, -np.inf], 0)  # replace infinite values with zeros

    # Plot clustergram without ax argument
    g = sns.clustermap(df_clustergram, method='average', cmap="coolwarm", standard_scale=1, figsize=(10, 8))
    
    # Convert the figure to a BytesIO object for Streamlit compatibility
    buf = io.BytesIO()
    plt.savefig(buf, format="png")
    buf.seek(0)
    return buf


def enrichr_analysis(gene_list, dataset):
    ENRICHR_URL = "http://amp.pharm.mssm.edu/Enrichr/addList"
    genes_str = "\n".join(gene_list)
    payload = {"list": (None, genes_str)}

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception("Error analyzing gene list")

    data = response.json()

    ENRICHR_URL = f'http://amp.pharm.mssm.edu/Enrichr/enrich?userListId={data["userListId"]}&backgroundType={dataset}'
    response = requests.get(ENRICHR_URL)
    if not response.ok:
        raise Exception("Error fetching enrichment results")

    num_genes = len(gene_list)  # Get the number of genes in the gene_list
    data = response.json()[dataset][0:num_genes]  # Dynamic slicing based on number of genes
    terms = [x[1] for x in data]
    scores = [-np.log10(x[2]) for x in data]

    return terms, scores

# # Function to perform enrichR analysis
# def enrichr_analysis(gene_list, dataset):
#     ENRICHR_URL = "http://amp.pharm.mssm.edu/Enrichr/addList"
#     genes_str = "\n".join(gene_list)
#     payload = {"list": (None, genes_str)}

#     response = requests.post(ENRICHR_URL, files=payload)
#     if not response.ok:
#         raise Exception("Error analyzing gene list")

#     data = response.json()

#     ENRICHR_URL = f'http://amp.pharm.mssm.edu/Enrichr/enrich?userListId={data["userListId"]}&backgroundType={dataset}'
#     response = requests.get(ENRICHR_URL)
#     if not response.ok:
#         raise Exception("Error fetching enrichment results")

#     data = response.json()[dataset][0:10]  # Getting the top 10 enriched terms
#     terms = [x[1] for x in data]
#     scores = [-np.log10(x[2]) for x in data]

#     return terms, scores


def plot_enrichr_results(terms, scores, title="Enrichr Results", widget_key=None):

    color_options = [
        "magma",
        "inferno",
        "plasma",
        "viridis",
        "cividis",
        "twilight",
        "twilight_shifted",
        "turbo",
        "Blues",
        "BrBG",
        "BuGn",
        "BuPu",
        "CMRmap",
        "GnBu",
        "Greens",
        "Greys",
        "OrRd",
        "Oranges",
        "PRGn",
        "PiYG",
        "PuBu",
        "PuBuGn",
        "PuOr",
        "PuRd",
        "Purples",
        "RdBu",
        "RdGy",
        "RdPu",
        "RdYlBu",
        "RdYlGn",
        "Reds",
        "Spectral",
        "Wistia",
        "YlGn",
        "YlGnBu",
        "YlOrBr",
        "YlOrRd",
        "afmhot",
        "autumn",
        "binary",
        "bone",
        "brg",
        "bwr",
        "cool",
        "coolwarm",
        "copper",
        "cubehelix",
        "flag",
        "gist_earth",
        "gist_gray",
        "gist_heat",
        "gist_ncar",
        "gist_rainbow",
        "gist_stern",
        "gist_yarg",
        "gnuplot",
        "gnuplot2",
        "gray",
        "hot",
        "hsv",
        "jet",
        "nipy_spectral",
        "ocean",
        "pink",
        "prism",
        "rainbow",
        "seismic",
        "spring",
        "summer",
        "terrain",
        "winter",
        "Accent",
        "Dark2",
        "Paired",
        "Pastel1",
        "Pastel2",
        "Set1",
        "Set2",
        "Set3",
        "tab10",
        "tab20",
        "tab20b",
        "tab20c",
        "magma_r",
        "inferno_r",
        "plasma_r",
        "viridis_r",
        "cividis_r",
        "twilight_r",
        "twilight_shifted_r",
        "turbo_r",
        "Blues_r",
        "BrBG_r",
        "BuGn_r",
        "BuPu_r",
        "CMRmap_r",
        "GnBu_r",
        "Greens_r",
        "Greys_r",
        "OrRd_r",
        "Oranges_r",
        "PRGn_r",
        "PiYG_r",
        "PuBu_r",
        "PuBuGn_r",
        "PuOr_r",
        "PuRd_r",
        "Purples_r",
        "RdBu_r",
        "RdGy_r",
        "RdPu_r",
        "RdYlBu_r",
        "RdYlGn_r",
        "Reds_r",
        "Spectral_r",
        "Wistia_r",
        "YlGn_r",
        "YlGnBu_r",
        "YlOrBr_r",
        "YlOrRd_r",
        "afmhot_r",
        "autumn_r",
        "binary_r",
        "bone_r",
        "brg_r",
        "bwr_r",
        "cool_r",
        "coolwarm_r",
        "copper_r",
        "cubehelix_r",
        "flag_r",
        "gist_earth_r",
        "gist_gray_r",
        "gist_heat_r",
        "gist_ncar_r",
        "gist_rainbow_r",
        "gist_stern_r",
        "gist_yarg_r",
        "gnuplot_r",
        "gnuplot2_r",
        "gray_r",
        "hot_r",
        "hsv_r",
        "jet_r",
        "nipy_spectral_r",
        "ocean_r",
        "pink_r",
        "prism_r",
        "rainbow_r",
        "seismic_r",
        "spring_r",
        "summer_r",
        "terrain_r",
        "winter_r",
        "Accent_r",
        "Dark2_r",
        "Paired_r",
        "Pastel1_r",
        "Pastel2_r",
        "Set1_r",
        "Set2_r",
        "Set3_r",
        "tab10_r",
        "tab20_r",
        "tab20b_r",
        "tab20c_r",
    ]
    selected_color = st.sidebar.selectbox(
        "Select Color Scheme", color_options, key=widget_key
    )

    plt.figure(
        figsize=(12, 6)
    )  # Increase the figure size to provide more space for labels
    sns.barplot(x=scores, y=terms, palette=selected_color)

    plt.title(title)
    plt.xlabel("-Log10(p-value)")
    plt.tight_layout()  # Ensuring labels are well displayed

    buf = BytesIO()
    plt.savefig(buf, format="tiff")
    buf.seek(0)
    plt.close()
    return buf


def modified_volcano_plot(
    df,
    x_col,
    y_col,
    gene_col,
    title,
    x_threshold_up,
    y_threshold_up,
    x_threshold_down,
    y_threshold_down,
    color_up,
    color_down,
):
    plt.figure(figsize=(10, 8))
    sns.scatterplot(data=df, x=x_col, y=y_col, color="grey", s=50)

    # Filter data for annotations and plots
    upregulated = df[(df[x_col] > x_threshold_up) & (df[y_col] > y_threshold_up)]
    downregulated = df[(df[x_col] < x_threshold_down) & (df[y_col] > y_threshold_down)]

    # Plotting and annotating upregulated genes
    sns.scatterplot(data=upregulated, x=x_col, y=y_col, color=color_up, s=50)
    for _, row in upregulated.iterrows():
        plt.annotate(row[gene_col], (row[x_col], row[y_col]), fontsize=9, alpha=0.7)

    # Plotting and annotating downregulated genes
    sns.scatterplot(data=downregulated, x=x_col, y=y_col, color=color_down, s=50)
    for _, row in downregulated.iterrows():
        plt.annotate(row[gene_col], (row[x_col], row[y_col]), fontsize=9, alpha=0.7)

    plt.axvline(x=x_threshold_up, color=color_up, linestyle="--")
    plt.axvline(x=x_threshold_down, color=color_down, linestyle="--")
    plt.axhline(y=y_threshold_up, color=color_up, linestyle="--")
    plt.axhline(y=y_threshold_down, color=color_down, linestyle="--")
    plt.title(title)

    buf = BytesIO()
    plt.savefig(buf, format="tiff")
    buf.seek(0)
    plt.close()
    return buf


def draw_venn(set1, set2, title="Venn Diagram"):
    plt.figure(figsize=(8, 8))
    venn2([set1, set2], set_labels=("MassSpec", "RNA-seq"))
    plt.title(title)
    buf = BytesIO()
    plt.savefig(buf, format="tiff")
    buf.seek(0)
    plt.close()
    return buf


def save_and_get_image_link(fig, filename):
    # Convert to Base64
    data = base64.b64encode(fig.getvalue()).decode()
    href = f'<a href="data:image/tiff;base64,{data}" download="{filename}">Download {filename}</a>'
    return href


def show_gene_lists(
    massspec_up, rna_up, massspec_down, rna_down, massspec_df=None, rna_df=None
):
    st.sidebar.title("Gene Lists")
    choice = st.sidebar.radio(
        "Choose a gene list",
        [
            "Upregulated in MassSpec",
            "Upregulated in RNA-seq",
            "Downregulated in MassSpec",
            "Downregulated in RNA-seq",
            "Common Upregulated Genes",
            "Common Downregulated Genes",
            "Spearman Correlation for Common Upregulated Genes",
            "Spearman Correlation for Common Downregulated Genes",
        ],
    )

    if choice in [
        "Upregulated in MassSpec",
        "Upregulated in RNA-seq",
        "Downregulated in MassSpec",
        "Downregulated in RNA-seq",
        "Common Upregulated Genes",
        "Common Downregulated Genes",
    ]:
        if choice == "Upregulated in MassSpec":
            data = list(massspec_up)
        elif choice == "Upregulated in RNA-seq":
            data = list(rna_up)
        elif choice == "Downregulated in MassSpec":
            data = list(massspec_down)
        elif choice == "Downregulated in RNA-seq":
            data = list(rna_down)
        elif choice == "Common Upregulated Genes":
            data = list(massspec_up.intersection(rna_up))
        elif choice == "Common Downregulated Genes":
            data = list(massspec_down.intersection(rna_down))

        # Displaying genes as a DataFrame
        df = pd.DataFrame(data, columns=["Gene Names"])
        st.table(df)

    elif choice == "Spearman Correlation for Common Upregulated Genes":
        up_common = massspec_up.intersection(rna_up)
        if massspec_df is not None and rna_df is not None and up_common:
            spearman_corr_up = massspec_df[massspec_df["Gene names"].isin(up_common)][
                "Log2 Ratio"
            ].corr(
                rna_df[rna_df["Gene names"].isin(up_common)]["Log2FC"],
                method="spearman",
            )
            st.write(
                f"Spearman correlation for common upregulated genes: {spearman_corr_up:.3f}"
            )
        else:
            st.warning("Datasets not provided or no common genes!")

    elif choice == "Spearman Correlation for Common Downregulated Genes":
        down_common = massspec_down.intersection(rna_down)
        if massspec_df is not None and rna_df is not None and down_common:
            spearman_corr_down = massspec_df[
                massspec_df["Gene names"].isin(down_common)
            ]["Log2 Ratio"].corr(
                rna_df[rna_df["Gene names"].isin(down_common)]["Log2FC"],
                method="spearman",
            )
            st.write(
                f"Spearman correlation for common downregulated genes: {spearman_corr_down:.3f}"
            )
        else:
            st.warning("Datasets not provided or no common genes!")


def get_download_link(buffer, filename, text):
    """Generate a link to download the image"""
    b64 = base64.b64encode(buffer.getvalue()).decode()
    return f'<a href="data:image/tiff;base64,{b64}" download="{filename}">{text}</a>'





def main():
    st.sidebar.title("Navigation")
    page_selection = st.sidebar.radio(
        "Choose a Page", ["Volcano and Venn Plots", "Gene Lists"]
    )

    # Initialization
    up_common_genes = set()
    down_common_genes = set()

    if page_selection == "Volcano and Venn Plots":
        st.title("Figures - Volcano and Venn Plots")
        st.audio("mp3.mp3", format="audio/mp3", start_time=0)

        massspec_file = st.file_uploader(
            "Drop MassSpec Dataset (.xlsx)", type=["xlsx"], key="massspec"
        )
        rna_file = st.file_uploader(
            "Drop RNA seq Dataset (.xlsx)", type=["xlsx"], key="rna"
        )

        massspec_upregulated_genes = set()
        rna_upregulated_genes = set()
        massspec_downregulated_genes = set()
        rna_downregulated_genes = set()

        if massspec_file and rna_file:
            massspec_df = pd.read_excel(massspec_file)
            rna_df = pd.read_excel(rna_file)

            # MassSpec volcano plot settings
            x_threshold_MS_up = st.slider(
                "Threshold for MassSpec Log2 Ratio (Upregulated)", 0.0, 5.0, 2.0
            )
            y_threshold_MS_up = st.slider(
                "Threshold for MassSpec -Log p-value (Upregulated)", 0.0, 5.0, 1.3
            )
            x_threshold_MS_down = st.slider(
                "Threshold for MassSpec Log2 Ratio (Downregulated)", -5.0, 0.0, -2.0
            )
            y_threshold_MS_down = st.slider(
                "Threshold for MassSpec -Log p-value (Downregulated)", 0.0, 5.0, 1.3
            )
            color_up_MS = st.color_picker(
                "Choose a color for Upregulated thresholds in MassSpec Volcano Plot",
                "#FF0000",
            )
            color_down_MS = st.color_picker(
                "Choose a color for Downregulated thresholds in MassSpec Volcano Plot",
                "#0000FF",
            )
            title_MS = st.text_input(
                "Title for MassSpec Volcano Plot", "Volcano Plot - MassSpec Data"
            )
            fig_MS = modified_volcano_plot(
                massspec_df,
                "Log2 Ratio",
                "˗Log p-value",
                "Gene names",
                title_MS,
                x_threshold_MS_up,
                y_threshold_MS_up,
                x_threshold_MS_down,
                y_threshold_MS_down,
                color_up_MS,
                color_down_MS,
            )

            # RNA-seq volcano plot settings
            x_threshold_RNA_up = st.slider(
                "Threshold for RNA-seq Log2 FC (Upregulated)", 0.0, 5.0, 2.0
            )
            y_threshold_RNA_up = st.slider(
                "Threshold for RNA-seq -Log p-value (Upregulated)", 0.0, 5.0, 1.3
            )
            x_threshold_RNA_down = st.slider(
                "Threshold for RNA-seq Log2 FC (Downregulated)", -5.0, 0.0, -2.0
            )
            y_threshold_RNA_down = st.slider(
                "Threshold for RNA-seq -Log p-value (Downregulated)", 0.0, 5.0, 1.3
            )
            color_up_RNA = st.color_picker(
                "Choose a color for Upregulated thresholds in RNA-seq Volcano Plot",
                "#FF0000",
            )
            color_down_RNA = st.color_picker(
                "Choose a color for Downregulated thresholds in RNA-seq Volcano Plot",
                "#0000FF",
            )
            title_RNA = st.text_input(
                "Title for RNA-seq Volcano Plot", "Volcano Plot - RNA-seq Data"
            )
            fig_RNA = modified_volcano_plot(
                rna_df,
                "Log2FC",
                "˗Log p-value",
                "Gene names",
                title_RNA,
                x_threshold_RNA_up,
                y_threshold_RNA_up,
                x_threshold_RNA_down,
                y_threshold_RNA_down,
                color_up_RNA,
                color_down_RNA,
            )

            # Venn diagrams for upregulated genes
            massspec_upregulated_genes = set(
                massspec_df[massspec_df["Log2 Ratio"] > x_threshold_MS_up]["Gene names"]
            )
            rna_upregulated_genes = set(
                rna_df[rna_df["Log2FC"] > x_threshold_RNA_up]["Gene names"]
            )
            fig_upregulated = draw_venn(
                massspec_upregulated_genes,
                rna_upregulated_genes,
                "Venn Diagram - Upregulated Genes",
            )

            # Venn diagrams for downregulated genes
            massspec_downregulated_genes = set(
                massspec_df[massspec_df["Log2 Ratio"] < x_threshold_MS_down][
                    "Gene names"
                ]
            )
            rna_downregulated_genes = set(
                rna_df[rna_df["Log2FC"] < x_threshold_RNA_down]["Gene names"]
            )
            fig_downregulated = draw_venn(
                massspec_downregulated_genes,
                rna_downregulated_genes,
                "Venn Diagram - Downregulated Genes",
            )

            # Display the plots directly in the Streamlit app
            st.image(
                fig_MS.getvalue(),
                caption="MassSpec Volcano Plot",
                use_column_width=True,
            )
            st.image(
                fig_RNA.getvalue(),
                caption="RNA-seq Volcano Plot",
                use_column_width=True,
            )
            st.image(
                fig_upregulated.getvalue(),
                caption="Venn Diagram - Upregulated Genes",
                use_column_width=True,
            )
            st.image(
                fig_downregulated.getvalue(),
                caption="Venn Diagram - Downregulated Genes",
                use_column_width=True,
            )

            # Create download buttons for each figure
            st.markdown(
                get_download_link(
                    fig_MS,
                    "massspec_volcano_plot.tiff",
                    "Download MassSpec Volcano Plot",
                ),
                unsafe_allow_html=True,
            )
            st.markdown(
                get_download_link(
                    fig_RNA,
                    "rna_seq_volcano_plot.tiff",
                    "Download RNA-seq Volcano Plot",
                ),
                unsafe_allow_html=True,
            )
            st.markdown(
                get_download_link(
                    fig_upregulated,
                    "upregulated_genes_venn_diagram.tiff",
                    "Download Venn Diagram of Upregulated Genes",
                ),
                unsafe_allow_html=True,
            )
            st.markdown(
                get_download_link(
                    fig_downregulated,
                    "downregulated_genes_venn_diagram.tiff",
                    "Download Venn Diagram of Downregulated Genes",
                ),
                unsafe_allow_html=True,
            )

            # After displaying download links for plots
            st.header("Enrichment Analysis")

            # Common genes between massspec and RNA for upregulated and downregulated genes
            up_common_genes = massspec_upregulated_genes.intersection(
                rna_upregulated_genes
            )
            down_common_genes = massspec_downregulated_genes.intersection(
                rna_downregulated_genes
            )

            # Dropdown for Enrichment Dataset
            datasets = [
                "KEGG_2019_Human",
                "KEGG_2019_Mouse",
                "GO_Biological_Process_2018",
                "Reactome_2016",
            ]
            # ... you can add more datasets from Enrichr here
            selected_dataset = st.selectbox(
                "Select a dataset for enrichment analysis:", datasets
            )

            # Display Enrichr results for common upregulated genes
            if up_common_genes:
                terms_up, scores_up = enrichr_analysis(up_common_genes, selected_dataset)
                fig_enrichr_up = plot_enrichr_results(
                    terms_up,
                    scores_up,
                    title="Enrichr Results for Common Upregulated Genes",
                    widget_key="color_scheme_up",
                )
                st.image(
                    fig_enrichr_up.getvalue(),
                    caption="Enrichr Results for Common Upregulated Genes",
                    use_column_width=True,
                )
                st.markdown(
                    get_download_link(
                        fig_enrichr_up,
                        "enrichr_upregulated_results.tiff",
                        "Download Enrichr Results for Upregulated Genes",
                    ),
                    unsafe_allow_html=True,
                )
                 # Display Enrichr clustergram for common upregulated genes
                clustergram_up = create_clustergram(list(up_common_genes), terms_up, scores_up)
                st.image(
                    clustergram_up,
                    caption="Clustergram for Common Upregulated Genes",
                    use_column_width=True,
                )
                st.markdown(
                    get_download_link(
                        clustergram_up,
                        "clustergram_upregulated.tiff",
                        "Download Clustergram for Upregulated Genes",
                    ),
                    unsafe_allow_html=True,
                )

            # Display Enrichr results for common downregulated genes
            if down_common_genes:
                terms_down, scores_down = enrichr_analysis(down_common_genes, selected_dataset)
                fig_enrichr_down = plot_enrichr_results(
                    terms_down,
                    scores_down,
                    title="Enrichr Results for Common Downregulated Genes",
                    widget_key="color_scheme_down",
                )
                st.image(
                    fig_enrichr_down.getvalue(),
                    caption="Enrichr Results for Common Downregulated Genes",
                    use_column_width=True,
                )
                st.markdown(
                    get_download_link(
                        fig_enrichr_down,
                        "enrichr_downregulated_results.tiff",
                        "Download Enrichr Results for Downregulated Genes",
                    ),
                    unsafe_allow_html=True,
                )
                # Display Enrichr clustergram for common downregulated genes
                clustergram_down = create_clustergram(list(down_common_genes), terms_down, scores_down)
                st.image(
                    clustergram_down,
                    caption="Clustergram for Common Downregulated Genes",
                    use_column_width=True,
                )
                st.markdown(
                    get_download_link(
                        clustergram_down,
                        "clustergram_downregulated.tiff",
                        "Download Clustergram for Downregulated Genes",
                    ),
                    unsafe_allow_html=True,
                )                
           # Check if session_state keys exist. If not, initialize them.
            if not hasattr(st.session_state, 'massspec_up'):
                st.session_state.massspec_up = set()

            if not hasattr(st.session_state, 'rna_up'):
                st.session_state.rna_up = set()

            if not hasattr(st.session_state, 'massspec_down'):
                st.session_state.massspec_down = set()

            if not hasattr(st.session_state, 'rna_down'):
                st.session_state.rna_down = set()

            if not hasattr(st.session_state, 'massspec_df'):
                st.session_state.massspec_df = pd.DataFrame()  # Initialize with an empty dataframe or any default value

            if not hasattr(st.session_state, 'rna_df'):
                st.session_state.rna_df = pd.DataFrame()  # Initialize with an empty dataframe or any default value


        # Storing gene sets in session state for use in the "Gene Lists" page
        st.session_state.massspec_up = massspec_upregulated_genes
        st.session_state.rna_up = rna_upregulated_genes
        st.session_state.massspec_down = massspec_downregulated_genes
        st.session_state.rna_down = rna_downregulated_genes


    elif page_selection == "Gene Lists":
        st.title("Gene Lists")

        # Check if the datasets are loaded and processed
        if ("massspec_df" in st.session_state) and ("rna_df" in st.session_state):
            # Display the gene lists
            st.subheader("MassSpec Upregulated Genes")
            st.table(st.session_state.massspec_up)

            st.subheader("RNA Upregulated Genes")
            st.table(st.session_state.rna_up)

            st.subheader("MassSpec Downregulated Genes")
            st.table(st.session_state.massspec_down)

            st.subheader("RNA Downregulated Genes")
            st.table(st.session_state.rna_down)

            # Common genes between massspec and RNA for upregulated and downregulated genes
            up_common_genes = st.session_state.massspec_up.intersection(st.session_state.rna_up)
            down_common_genes = st.session_state.massspec_down.intersection(st.session_state.rna_down)

            # Displaying common upregulated genes
            if up_common_genes:
                st.subheader("Common Upregulated Genes")
                st.table(list(up_common_genes))

            # Displaying common downregulated genes
            if down_common_genes:
                st.subheader("Common Downregulated Genes")
                st.table(list(down_common_genes))

# #             # Spearman correlation for upregulated genes
#             massspec_subset_up = st.session_state.masspec_df[st.session_state.masspec_df["Gene names"].isin(up_common_genes)][["Gene names", "Log2 Ratio"]].set_index("Gene names")
#             rna_subset_up = st.session_state.rna_df[st.session_state.rna_df["Gene names"].isin(up_common_genes)][["Gene names", "Log2FC"]].set_index("Gene names")

#             # Merge datasets for common genes, ensuring same order
#             merged_up = massspec_subset_up.join(rna_subset_up, how='inner')
#             spearman_corr_up = merged_up["Log2 Ratio"].corr(merged_up["Log2FC"], method='spearman')

#             st.subheader("Spearman Correlation for Common Upregulated Genes")
#             st.write(spearman_corr_up)


    else:
        st.warning("Please upload datasets in the 'Volcano and Venn Plots' section first.")



# Run the Streamlit app
if __name__ == '__main__':
    main()

