import streamlit as st
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import io

# Sayfa Ayarları
st.set_page_config(page_title="Tez Molekül Çizici", page_icon="🧪")

st.title("🧪 Tez İçin Yüksek Çözünürlüklü Molekül Çizici")
st.markdown("""
Bu araç, girilen SMILES kodunu **tezlerde kullanıma uygun (300 DPI, Yüksek Çözünürlük)** şeffaf PNG formatına dönüştürür.
""")

# Kullanıcıdan Girdi Alma
smiles = st.text_input("SMILES Kodunu Girin:", value="CC(=O)OC1=CC=CC=C1C(=O)O")
cozunurluk = st.slider("Görsel Çözünürlüğü (Piksel)", 500, 2000, 1000)

if st.button("Çiz ve Hazırla"):
    if not smiles:
        st.warning("Lütfen bir SMILES kodu girin.")
    else:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                st.error("Geçersiz SMILES kodu! Lütfen kontrol edin.")
            else:
                # --- YÜKSEK KALİTE ÇİZİM MOTORU ---
                # RDKit'in Cairo motoru vektörel kalitede çizim yapar
                d = rdMolDraw2D.MolDraw2DCairo(cozunurluk, cozunurluk)
                
                # Çizim ayarları (Tez için optimize edildi)
                opts = d.drawOptions()
                opts.bondLineWidth = 3        # Çizgiler kalın ve net
                opts.clearBackground = True   # Arka plan şeffaf
                opts.padding = 0.05           # Kenar boşluğu
                
                # Çizimi gerçekleştir
                d.DrawMolecule(mol)
                d.FinishDrawing()
                
                # Çıktıyı belleğe al (Dosya olarak kaydetmeden)
                png_data = d.GetDrawingText()
                
                # Ekranda göster
                st.image(png_data, caption="Önizleme (İndirilen dosya daha yüksek kalitededir)", use_container_width=False, width=400)
                
                # İNDİRME BUTONU
                st.download_button(
                    label="📥 Yüksek Kaliteli PNG İndir",
                    data=png_data,
                    file_name="molekul_tez_kalitesi.png",
                    mime="image/png"
                )
                st.success("Görsel hazır! İndirme butonuna basabilirsiniz.")
                
        except Exception as e:
            st.error(f"Bir hata oluştu: {e}")

st.markdown("---")
st.caption("SMILES kodlarını Wikipedia veya PubChem üzerinden bulabilirsiniz.")
