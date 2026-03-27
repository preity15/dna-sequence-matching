import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# Page Config + Styling
# -------------------------------
st.set_page_config(page_title="DNA Matcher", page_icon="🧬", layout="centered")

st.markdown("""
    <style>
    .stApp {
        background: linear-gradient(to right, #0f2027, #203a43, #2c5364);
        color: white;
    }
    h1, h2, h3 {
        text-align: center;
    }
    </style>
""", unsafe_allow_html=True)

# -------------------------------
# Header + Static DNA Image
# -------------------------------
st.markdown("<h1>🧬 DNA Sequence Matching System</h1>", unsafe_allow_html=True)

st.image("https://upload.wikimedia.org/wikipedia/commons/8/87/DNA_chemical_structure.svg", width=200)

st.write("### Compare DNA sequences using Needleman-Wunsch Algorithm")

# -------------------------------
# Needleman-Wunsch Algorithm
# -------------------------------
def needleman_wunsch(seq1, seq2):
    n, m = len(seq1), len(seq2)

    match = 1
    mismatch = -1
    gap = -2

    dp = [[0]*(m+1) for _ in range(n+1)]

    for i in range(n+1):
        dp[i][0] = i * gap
    for j in range(m+1):
        dp[0][j] = j * gap

    for i in range(1, n+1):
        for j in range(1, m+1):

            if seq1[i-1] == seq2[j-1]:
                score = match
            else:
                score = mismatch

            dp[i][j] = max(
                dp[i-1][j-1] + score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )

    # Traceback
    aligned1, aligned2 = "", ""
    i, j = n, m

    while i > 0 and j > 0:
        current = dp[i][j]
        diag = dp[i-1][j-1]
        up = dp[i-1][j]
        left = dp[i][j-1]

        if seq1[i-1] == seq2[j-1]:
            score = match
        else:
            score = mismatch

        if current == diag + score:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        elif current == up + gap:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    while i > 0:
        aligned1 = seq1[i-1] + aligned1
        aligned2 = "-" + aligned2
        i -= 1

    while j > 0:
        aligned1 = "-" + aligned1
        aligned2 = seq2[j-1] + aligned2
        j -= 1

    return dp, dp[n][m], aligned1, aligned2


# -------------------------------
# Helper Functions
# -------------------------------
def show_alignment(a1, a2):
    line = ""
    for i in range(len(a1)):
        if a1[i] == a2[i]:
            line += "|"
        else:
            line += " "
    return line


def calculate_similarity(a1, a2):
    matches = sum(1 for i in range(len(a1)) if a1[i] == a2[i])
    return (matches / len(a1)) * 100


def is_valid_dna(seq):
    return all(c in "ATGC" for c in seq)


# -------------------------------
# Input Section
# -------------------------------
seq1 = st.text_input("Enter Sequence 1")
seq2 = st.text_input("Enter Sequence 2")

# -------------------------------
# Button Action
# -------------------------------
if st.button("🔍 Align Sequences"):

    if seq1 and seq2:
        seq1 = seq1.upper()
        seq2 = seq2.upper()

        if not is_valid_dna(seq1) or not is_valid_dna(seq2):
            st.error("❌ Use only A, T, G, C")
        else:
            dp, score, a1, a2 = needleman_wunsch(seq1, seq2)

            match_line = show_alignment(a1, a2)
            similarity = calculate_similarity(a1, a2)

            # -------------------------------
            # Output Section
            # -------------------------------
            st.subheader("📊 Alignment Score")
            st.success(score)

            st.subheader("🧬 Sequence Alignment")

            st.markdown(f"""
            <pre style="font-size:18px; color:#00FFAA;">
{a1}
{match_line}
{a2}
            </pre>
            """, unsafe_allow_html=True)

            st.subheader("📈 Similarity")
            st.info(f"{similarity:.2f}%")

            # -------------------------------
            # DP Matrix Visualization (IMPROVED)
            # -------------------------------
            st.subheader("📉 DP Matrix Visualization")

            fig, ax = plt.subplots()
            cax = ax.imshow(dp)

            # Show values inside cells
            for i in range(len(dp)):
                for j in range(len(dp[0])):
                    ax.text(j, i, dp[i][j], ha="center", va="center", fontsize=8)

            ax.set_title("Dynamic Programming Matrix")
            ax.set_xlabel("Sequence 2")
            ax.set_ylabel("Sequence 1")

            plt.colorbar(cax)

            st.pyplot(fig)

            # -------------------------------
            # Download Report
            # -------------------------------
            report = f"""
DNA Sequence Matching Report

Sequence 1: {seq1}
Sequence 2: {seq2}

Alignment Score: {score}
Similarity: {similarity:.2f}%

Alignment:
{a1}
{match_line}
{a2}
"""

            st.download_button(
                label="📄 Download Report",
                data=report,
                file_name="dna_report.txt",
                mime="text/plain"
            )

    else:
        st.warning("⚠️ Please enter both sequences")