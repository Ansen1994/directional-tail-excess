import os, struct, sys, urllib.request, zlib

URL = "https://zenodo.org/records/15799566/files/Source%20Data.zip?download=1"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "raw")
TAIL = 16 * 1024 * 1024

TARGETS = [
    "Source Data/Figure2/nmf_basis.pickle",
    "Source Data/Figure2/MP_Programs.gmt",
    "Source Data/Figure2/MP_anno_map.pickle",
    "Source Data/Figure2/MetaProgram_Annotation.csv",
    "Source Data/Figure1/GEX_OBS_Lineage.csv",
    "Source Data/Figure2/Annotation.csv",
    "Source Data/Figure2/signature_scores.csv",
    "Source Data/Figure2/MP_order_color.pickle",
    "Source Data/Figure1/gex_cnv_SCEVAN_karyotype.csv",
    "Source Data/Figure5/GEX_OBS_Cellstate.csv",
    "Source Data/Figure1/embeddings.pickle",
]

req = urllib.request.Request(URL, method="HEAD")
with urllib.request.urlopen(req, timeout=60) as r:
    final_url, size = r.geturl(), int(r.headers["Content-Length"])

def fetch(a, b):
    rq = urllib.request.Request(final_url, headers={"Range": f"bytes={a}-{b}"})
    with urllib.request.urlopen(rq, timeout=600) as r:
        return r.read()

start = max(0, size - TAIL)
tail = fetch(start, size - 1)

eocd = tail.rfind(b"PK\x05\x06")
n, cd_size = struct.unpack("<H", tail[eocd+10:eocd+12])[0], struct.unpack("<I", tail[eocd+12:eocd+16])[0]
cd_off = struct.unpack("<I", tail[eocd+16:eocd+20])[0]
if cd_off == 0xFFFFFFFF:
    loc = tail.rfind(b"PK\x06\x07")
    z64 = struct.unpack("<Q", tail[loc+8:loc+16])[0]
    rec = tail[z64-start: z64-start+56]
    n = struct.unpack("<Q", rec[32:40])[0]
    cd_size = struct.unpack("<Q", rec[40:48])[0]
    cd_off = struct.unpack("<Q", rec[48:56])[0]

cd = tail[cd_off-start: cd_off-start+cd_size]

entries = {}
i = 0
while i < len(cd) - 46:
    if cd[i:i+4] != b"PK\x01\x02":
        i += 1; continue
    flag, = struct.unpack("<H", cd[i+8:i+10])
    csize, usize = struct.unpack("<II", cd[i+20:i+28])
    nlen, elen, clen = struct.unpack("<HHH", cd[i+28:i+34])
    lho, = struct.unpack("<I", cd[i+42:i+46])
    name = cd[i+46:i+46+nlen].decode("utf-8", "replace")
    extra = cd[i+46+nlen:i+46+nlen+elen]
    if usize == 0xFFFFFFFF or csize == 0xFFFFFFFF or lho == 0xFFFFFFFF:
        j = 0
        while j + 4 <= len(extra):
            hid, hsz = struct.unpack("<HH", extra[j:j+4])
            if hid == 0x0001:
                b = extra[j+4:j+4+hsz]; k = 0
                if usize == 0xFFFFFFFF: usize = struct.unpack("<Q", b[k:k+8])[0]; k += 8
                if csize == 0xFFFFFFFF: csize = struct.unpack("<Q", b[k:k+8])[0]; k += 8
                if lho == 0xFFFFFFFF: lho = struct.unpack("<Q", b[k:k+8])[0]; k += 8
                break
            j += 4 + hsz
    entries[name] = (lho, csize, usize)
    i += 46 + nlen + elen + clen

missing = [t for t in TARGETS if t not in entries]
if missing:
    print("NOT FOUND:"); [print("  ", m) for m in missing]

os.makedirs(os.path.join(OUT, "nc2025"), exist_ok=True)
for t in TARGETS:
    if t not in entries:
        continue
    lho, csize, usize = entries[t]
    hdr = fetch(lho, lho + 29)
    nlen, elen = struct.unpack("<HH", hdr[26:30])
    data_start = lho + 30 + nlen + elen
    blob = fetch(data_start, data_start + csize - 1)
    try:
        raw = zlib.decompressobj(-15).decompress(blob)
    except Exception as e:
        print(f"decompress FAILED {t}: {e}"); continue
    fn = os.path.join(OUT, "nc2025", os.path.basename(t))
    if os.path.exists(fn) and os.path.getsize(fn) == usize:
        print(f"SKIP {t}"); continue
    with open(fn, "wb") as f:
        f.write(raw)
    print(f"OK {usize/1e6:8.2f} MB -> {os.path.basename(t)}")
