import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import SQLAlchemyError

# Ruta al archivo CSV y detalles de la base de datos MySQL
csv_file = 'heatmap.csv'
mysql_username = 'root'
mysql_password = 'Cim12345!'
mysql_host = 'localhost'
mysql_dbname = 'datos'

# Crear una conexión a la base de datos MySQL
engine = create_engine(f'mysql+mysqlconnector://{mysql_username}:{mysql_password}@{mysql_host}/{mysql_dbname}')
Session = sessionmaker(bind=engine)
session = Session()

try:
    # SENSOR VALUES
    dtypes1 = {'temp': 'float', 'mltss': 'int', 'sludge_prod': 'float'}
    df1 = pd.read_csv(csv_file, skiprows=1, dtype=dtypes1, names=['temp', 'mltss', 'sludge_prod'])
    df1.to_sql('heatmap', con=engine, if_exists='append', index=False)
    print(f'Datos importados exitosamente a la tabla "sensor_datos" en MySQL desde el archivo CSV.')

    # Commit the transaction
    session.commit()
except SQLAlchemyError as e:
    # Rollback the transaction on error
    session.rollback()
    print(f'Error al importar datos: {str(e)}')
finally:
    # Close the session
    session.close()
    print("Sesión cerrada")
