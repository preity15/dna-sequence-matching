# 🧬 DNA Sequence Matching System

A bioinformatics-based web application that performs DNA sequence alignment using the **Needleman-Wunsch Algorithm**. This project demonstrates the application of **Dynamic Programming** in solving real-world problems like sequence matching and similarity analysis.

## Features

- 🔬 Global DNA sequence alignment
- 📊 Alignment score calculation
- 📈 Similarity percentage computation
- 🧬 Visual alignment with match indicators
- 📉 Dynamic Programming (DP) matrix visualization
- 📄 Downloadable alignment report
- 🎨 Clean and interactive UI using Streamlit


## Tech Stack

- **Language:** Python  
- **Frontend/UI:** Streamlit  
- **Libraries:** NumPy, Matplotlib  
- **Concepts Used:**  
  - Dynamic Programming  
  - String Matching  
  - Optimization Algorithms  

## How It Works

1. User inputs two DNA sequences (A, T, G, C)
2. The system applies the **Needleman-Wunsch algorithm**
3. A DP matrix is constructed to compute optimal alignment
4. Traceback is performed to generate aligned sequences
5. Results are displayed with:
   - Alignment score  
   - Matching visualization  
   - Similarity percentage  


## ▶ Run Locally

pip install -r requirements.txt
streamlit run app.py