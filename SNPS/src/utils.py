import json
import os

def guardar_json(datos, ruta):
    os.makedirs(os.path.dirname(ruta), exist_ok=True)

    with open(ruta, "w") as f:
        json.dump(datos, f, indent=4)