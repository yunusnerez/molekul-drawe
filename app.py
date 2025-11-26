import streamlit as st
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import requests

# Sayfa Ayarları
st.set_page_config(page_title="Molekül Çizici", page_icon="🧪")

gizleme_kodu = """
<style>
    /* Sağ üstteki 'Deploy', 'Fork', 'Settings' menüsünü gizle */
    [data-testid="stToolbar"] {
        visibility: hidden;
        display: none;
    }
    
    /* Üstteki renkli şeridi ve header alanını gizle */
    .stApp > header {
        visibility: hidden;
        display: none;
    }
    
    /* En alttaki 'Made with Streamlit' yazısını gizle */
    footer {
        visibility: hidden;
        display: none;
    }
    
    /* Sayfanın üst kısmındaki boşluğu biraz azalt (header gidince boşluk kalmasın) */
    .block-container {
        padding-top: 2rem;
    }
</style>
"""
st.markdown(gizleme_kodu, unsafe_allow_html=True)

# --- FONKSİYON: PubChem'den Veri Çekme ---
def get_smiles_from_name(molecule_name):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"
    url = f"{base_url}/{molecule_name}/property/IsomericSMILES/TXT"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            return response.text.strip()
        else:
            return None
    except:
        return None

# --- ARAYÜZ ---
st.title("🧪 Akıllı Molekül Çizici (PubChem Entegreli)")
st.markdown("""
Molekülün **İngilizce adını** yazın (örn: *Aspirin, Ibuprofen, Caffeine*), sistem otomatik olarak SMILES kodunu bulup kutuya yazacaktır.
""")

# Oturum Durumu Başlatma (Hafıza)
if "smiles_entry" not in st.session_state:
    st.session_state.smiles_entry = "CC(=O)OC1=CC=CC=C1C(=O)O" # Varsayılan: Aspirin
if "search_result" not in st.session_state:
    st.session_state.search_result = None

# 1. ARAMA BÖLÜMÜ
col_search1, col_search2 = st.columns([3, 1])
with col_search1:
    search_name = st.text_input("Molekül Adı ile Ara (İngilizce):", placeholder="Örn: Cholesterol")
with col_search2:
    st.write("")
    st.write("") 
    if st.button("🔍 Bul"):
        if search_name:
            with st.spinner("Aranıyor..."):
                found_smiles = get_smiles_from_name(search_name)
                if found_smiles:
                    st.session_state.search_result = found_smiles
                    st.success(f"Bulundu: {search_name}")
                else:
                    st.session_state.search_result = None
                    st.error("Bulunamadı! İsmi İngilizce yazdığınızdan emin olun.")
        else:
            st.warning("Lütfen bir isim yazın.")

# Arama Sonucu ve Seçme Alanı
if st.session_state.search_result:
    st.info(f"Bulunan SMILES: {st.session_state.search_result}")
    if st.button("✅ Bu Molekülü Seç"):
        st.session_state.smiles_entry = st.session_state.search_result
        st.session_state.search_result = None # Seçtikten sonra sonucu temizle
        st.rerun()

st.markdown("---")

# 2. ÇİZİM BÖLÜMÜ
# key="smiles_entry" diyerek bu kutuyu yukarıdaki session_state'e bağladık.
# Artık yukarıda st.session_state.smiles_entry değişince burası da değişecek.
smiles = st.text_input("SMILES Kodu (Otomatik Dolatır veya Düzenleyebilirsiniz):", key="smiles_entry")

cozunurluk = st.slider("Görsel Çözünürlüğü (Piksel)", 500, 2000, 1000)

if st.button("🎨 Çizimi Oluştur"):
    if not smiles:
        st.warning("Lütfen SMILES kodu girin.")
    else:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                st.error("Geçersiz SMILES kodu!")
            else:
                # Çizim Motoru
                d = rdMolDraw2D.MolDraw2DCairo(cozunurluk, cozunurluk)
                opts = d.drawOptions()
                opts.bondLineWidth = 3
                opts.clearBackground = True
                opts.padding = 0.05
                
                d.DrawMolecule(mol)
                d.FinishDrawing()
                
                png_data = d.GetDrawingText()
                
                col1, col2 = st.columns([1, 1])
                with col1:
                    st.image(png_data, caption="Önizleme", use_container_width=False, width=400)
                with col2:
                    st.write("### Hazır! 👇")
                    st.download_button(
                        label="📥 PNG Olarak İndir",
                        data=png_data,
                        file_name="molekul_tez.png",
                        mime="image/png"
                    )
        except Exception as e:
            st.error(f"Hata: {e}")
