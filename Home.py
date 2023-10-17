import streamlit as st
from PIL import Image
import io

# Open the image file
with open("Logo.png", "rb") as image_file:
    image = Image.open(io.BytesIO(image_file.read()))

# Fun Colors
primary_color = "#ff6b6b"  # A bright, playful red
secondary_color = "#48dbfb"  # A bright, playful blue

# Custom CSS for a Fun Design
custom_css = f"""
    <style>
        /* Light background */
        body {{
            background-color: #f9f9f9;
        }}
        /* Fun Header Styling */
        .stMarkdown h1 {{
            font-family: 'Arial Rounded MT Bold', sans-serif;
            border-bottom: 3px solid {primary_color};
            padding: 10px 16px;
            border-radius: 15px;
            background: #ffffff;
            box-shadow: 3px 3px 10px lightgrey;
        }}
        /* Fun Subheader Styling */
        .stMarkdown h2 {{
            font-family: 'Arial Rounded MT Bold', sans-serif;
            color: {secondary_color};
            padding: 5px;
            background: #ffffff;
            border-radius: 10px;
            box-shadow: 2px 2px 5px lightgrey;
        }}
        /* Links with hover effect */
        .stMarkdown a {{
            color: {primary_color};
            text-decoration: none;
            transition: color .3s ease-in-out;
        }}
        .stMarkdown a:hover {{
            color: {secondary_color};
        }}
        /* Images styling */
        .stImage img {{
            border: 2px solid {primary_color};
            border-radius: 20px;
            box-shadow: 4px 4px 10px lightgrey;
        }}
        /* Fun Sidebar styling */
        .sidebar .sidebar-content {{
            background-color: #f9f9f9;
            border: none;
            box-shadow: inset 2px 2px 8px lightgrey;
        }}
    </style>
"""

st.markdown(custom_css, unsafe_allow_html=True)

# Header with Emoji and Title
st.markdown(f"<h1 style='color: {primary_color};'>🧬 Python tool for Proteo-Transcriptomic Analysis</h1>", unsafe_allow_html=True)
st.write("🌐 Welcome to the Proteo-Transcriptomic Streamlit web tool. Navigate through the sections to explore the functionalities!")

# Introduction
st.markdown(f"<h2>📘 Introduction</h2>", unsafe_allow_html=True)
st.write("""
Discover the unique functionalities of the Proteo-transcriptomics webtool kit, designed to unearth correlations between Protein expression datasets (Mass-Spectrometry reads) and RNA sequencing datasets (RNA-seq reads). In the field of molecular biology and genomics, a web tool for correlating proteomic and transcriptomic data is a valuable asset. It facilitates deeper insights into gene expression, identification of regulatory mechanisms, and biomarker discovery. This efficient automation accelerates research, enhances credibility, and expedites scientific discoveries, bridging the gap between genomic information and protein functionality for researchers
""")

# Image with shadow
st.image("Diag.png", use_column_width=True, caption="Sample visualization")

# File Prerequisites
st.markdown(f"<h2>📁 File pre-requisites</h2>", unsafe_allow_html=True)
st.markdown("""
- **MassSpec dataset**: Excel file (.xlsx) with 'Gene names', 'Log2 Ratio', '-Log p-value'.
- **RNA-Seq dataset**: Excel file (.xlsx) with 'Gene names', 'Log2FC', '-Log p-value'.
""", unsafe_allow_html=True)

# Key Features
st.markdown(f"<h2>🚀 Key Features</h2>", unsafe_allow_html=True)
st.write("""
Delve deep into the features offered:
1. 📊 Customizable volcano plots for both DEP and DEG datasets.
2. 📝 Lists of regulated proteins and transcript gene names.
3. 🤝 Venn Diagram of overlapping regulated genes.
4. 📈 Enrichr Analysis of common upregulated and downregulated genes.
""")

# Support Section
st.markdown(f"<h2>📞 Support</h2>", unsafe_allow_html=True)
st.write("""
For queries, reach out to:
- 📩 Look Zhuojian: [zhuojian.look@u.duke.nus.edu](mailto:zhuojian.look@u.duke.nus.edu)
- 📩 Sidra Mohamed Yaqoob: [sidra.yaqoob@u.duke.nus.edu](mailto:sidra.yaqoob@u.duke.nus.edu)
- 📩 Yadanar Than Naing: [e0869002@u.duke.nus.edu](mailto:e0869002@u.duke.nus.edu)
""")

# Sidebar
st.sidebar.markdown("<div style='height:450px;'></div>", unsafe_allow_html=True)
st.sidebar.image(image, use_column_width=True, caption="Official Logo")
