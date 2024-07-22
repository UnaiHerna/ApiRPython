import rpy2.robjects as robjects
import uvicorn
from fastapi import FastAPI, Depends, HTTPException
from sqlalchemy import select
from sqlalchemy.orm import Session
from starlette.middleware.cors import CORSMiddleware
from starlette.requests import Request
from starlette.responses import HTMLResponse
from starlette.templating import Jinja2Templates

from db.connector import get_db
from db.models import Heatmap
from db.redis_client import get_cached_response, set_cached_response

app = FastAPI()

# Configura CORS
app.add_middleware(
    # Permite todos los orígenes. Cambia esto para producción.
    # Permite todos los métodos (GET, POST, etc.)
    # Permite todos los encabezadosm
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


templates = Jinja2Templates(directory="web")


@app.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.get("/heatmap", response_class=HTMLResponse)
async def read_heatmap(request: Request):
    return templates.TemplateResponse("heatmap.html", {"request": request})


@app.get("/r/steadystate")
async def steadystate(mltss_sp: float, so_aer_sp: float, q_int: float, tss_eff_sp: float, temp: float):
    cache_key = f"datos_steady_{mltss_sp, so_aer_sp, q_int, tss_eff_sp, temp}"
    cached_data = get_cached_response(cache_key)
    if cached_data:
        return cached_data
    try:
        #inicio = time.time()

        r = robjects.r
        r['source']('f_steady_state_wwtp.R')

        # Cargamos la función que hemos definido en R
        function_r = robjects.globalenv['f_steady_state_wwtp']

        # Invocamos la función de R y obtenemos el resultado
        resultado = function_r(mltss_sp, so_aer_sp, q_int, tss_eff_sp, temp)

        result_dict = {}

        # Manejar el resultado: si es una lista (ListVector) o una cadena (Str)
        if isinstance(resultado, robjects.vectors.ListVector):
            for i in range(len(resultado)):
                if resultado.names[i] == "v_conc_anx" or resultado.names[i] == "v_conc_aer":
                    nested_dict = {}
                    for j in range(len(resultado[i])):
                        nested_dict[resultado[i].names[j]] = resultado[i][j]
                    result_dict[resultado.names[i]] = nested_dict
                else:
                    if len(resultado[i]) < 2:
                        result_dict[resultado.names[i]] = float(resultado[i][0])
                    else:
                        result_dict[resultado.names[i]] = list(resultado[i])
        elif isinstance(resultado, robjects.vectors.StrVector):
            result_dict["result"] = list(resultado)
        else:
            raise ValueError("El tipo de resultado devuelto por la función R no es compatible.")

        #fin = time.time()

        #print("Tiempo de ejecución: ", fin - inicio, "segundos")

        set_cached_response(cache_key, result_dict)
        return result_dict

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/datos/heatmap")
async def heatmap(db: Session = Depends(get_db)):
    cache_key = f"datos_heatmap_{db}"
    cached_data = get_cached_response(cache_key)
    if cached_data:
        return cached_data
    query = select(Heatmap).order_by(Heatmap.id)
    datos = db.execute(query).scalars().all()

    datos_dict = [dato.to_dict() for dato in datos]

    set_cached_response(cache_key, datos_dict)
    return datos_dict
