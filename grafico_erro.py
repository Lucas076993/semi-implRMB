import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


df = pd.read_csv("error.csv")


colunas = df.columns.tolist()
colunas.remove("t")

plt.figure(figsize=(12,6))
plt.title("Medidas de erro em relação ao tempo")
plt.xlabel("Tempo (t)")
plt.ylabel("Erro percentual (%)")

for i in colunas:  
  plt.scatter(df["t"], df[i], label=i, s=1)

plt.legend()
plt.savefig("erros.png")
