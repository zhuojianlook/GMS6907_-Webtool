import streamlit as st

def app():
    st.title("Thank You!")
    
    with open("fireworks.html", "r") as f:
        st.components.v1.html(f.read(), height=480, width=720)

if __name__ == '__main__':
    app()
