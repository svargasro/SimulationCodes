#ESTE ES PARA DIRECCIÓN

import pandas as pd
import numpy as np

localidades = np.array(('Movil_7ma','MinAmbiente','Puente_Aranda', 'Centro_de_Alto_Rendimiento'))

#localidades = np.array(('Centro_de_Alto_Rendimiento','3'))

for localidad in localidades:
    print(localidad)
    globaldirection = np.full(240, -1.0)


    for i in range(0,10):

      try:
        df = pd.read_excel(f'./datosDireccionViento/{localidad}/{localidad}_fromGraph_18-05-2024({i}).xlsx')
        print("try",i)
      except:
        df = pd.read_excel(f'./datosDireccionViento/{localidad}/{localidad}_fromGraph_18-05-2024.xlsx')
        print("except",i)

      dir = np.array(df)
      trans = np.transpose(dir)
      direction = trans[1][2:]*(np.pi/180)


      direction = np.array([np.nan if x == 'NaN' else x for x in direction], dtype=float)
      nan_indices = np.isnan(direction)

      for j in range(0, len(direction)):
          if nan_indices[j]:
              direction[j] = -1.0

      startIndex = i*(23+1)
      finalIndex= (i+1)*(23+1)
      print(len(direction))
      print(len(globaldirection[startIndex:finalIndex]))
      globaldirection[startIndex:finalIndex] = direction

    np.savetxt(f'./datosTXT/direccionViento_{localidad}.txt', globaldirection)



for localidad in localidades:
    print(localidad)
    globalVelocity = np.full(240, -1.0)

    for i in range(0,10):

      try:
        df = pd.read_excel(f'./datosVelocidadViento/{localidad}/{localidad}_fromGraph_18-05-2024({i}).xlsx')
        print("try",i)
      except:
        df = pd.read_excel(f'./datosVelocidadViento/{localidad}/{localidad}_fromGraph_18-05-2024.xlsx')
        print("except",i)

      vel = np.array(df)
      trans = np.transpose(vel)
      velocity = trans[1][2:]*(5*np.sqrt(3)/(5148))
      velocity = np.array([np.nan if x == 'NaN' else x for x in velocity], dtype=float)
      nan_indices = np.isnan(velocity)


      for j in range(0, len(velocity)):
          if nan_indices[j]:
              velocity[j] = -1

      startIndex = i*(23+1)
      finalIndex= (i+1)*(23+1)
      print(len (globalVelocity[startIndex:finalIndex]))
      globalVelocity[startIndex:finalIndex] = velocity

    np.savetxt(f'./datosTXT/velocidadViento_{localidad}.txt', globalVelocity)
