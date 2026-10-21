import re
import pandas as pd

# 1. Leggi il file
with open("Fucino_Lario_Sat_Access.txt", "r", encoding="utf-8") as f:
    text = f.read()

# 2. Dividi in due blocchi usando le intestazioni reali del file
parts = re.split(r"Satellite1 to Lario", text)
fucino_block = parts[0]   # Prima parte
lario_block = parts[1]    # Seconda parte

# 3. Regex per estrarre i dati
pattern = re.compile(
    r"(\d{1,2} \w+ \d{4})\s+(\d{2}:\d{2}:\d{2}\.\d+)\s+([0-9.]+)\s+([0-9.\-]+)\s+([0-9.]+)"
)

def parse_block(block):
    rows = []
    for match in pattern.finditer(block):
        date = match.group(1)
        time = match.group(2)
        az = float(match.group(3))
        el = float(match.group(4))
        rng = float(match.group(5))
        rows.append((date, time, az, el, rng))
    return rows

# 4. Estrai i dati
fucino_data = parse_block(fucino_block)
lario_data = parse_block(lario_block)

# 5. Converti in DataFrame
df_fucino = pd.DataFrame(fucino_data, columns=["Date","Time","Azimuth","Elevation","Range"])
df_lario = pd.DataFrame(lario_data, columns=["Date","Time","Azimuth","Elevation","Range"])

# 6. Esporta in CSV ordinati
df_fucino.to_csv("fucino_satellite.csv", index=False)
df_lario.to_csv("lario_satellite.csv", index=False)

print("CSV generati: fucino_satellite.csv e lario_satellite.csv")

