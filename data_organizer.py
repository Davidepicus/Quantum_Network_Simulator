import re
import pandas as pd

# 1. Leggi il file
with open("Lijiang-Micius_Delingha-Micius_Visibility_Windows.txt", "r", encoding="utf-8") as f:
    text = f.read()

# 2. Dividi in due blocchi
parts = re.split(r"Delingha to Micius", text)
lijiang_block = parts[0]
delingha_block = parts[1]

# 3. Regex per estrarre i dati
pattern = re.compile(
    r"(\d{1,2} \w+ \d{4})\s+(\d{2}:\d{2}:\d{2}\.\d+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)"
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
lijiang_data = parse_block(lijiang_block)
delingha_data = parse_block(delingha_block)

# 5. Converti in DataFrame
df_lijiang = pd.DataFrame(lijiang_data, columns=["Date","Time","Azimuth","Elevation","Range"])
df_delingha = pd.DataFrame(delingha_data, columns=["Date","Time","Azimuth","Elevation","Range"])

# 6. Esporta in CSV ordinati
df_lijiang.to_csv("lijiang_micius.csv", index=False)
df_delingha.to_csv("delingha_micius.csv", index=False)

print("CSV generati: lijiang_micius.csv e delingha_micius.csv")
