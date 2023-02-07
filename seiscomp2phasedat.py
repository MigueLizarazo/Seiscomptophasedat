'''
author: mlizarazo

La rutina seiscomp2phasedat.py que escribe los archivos phase.dat y station.dat (aptos como entrada de ph2dt (subrutina de hypoDD) y hypoDD, asi
como sus archivos analogos phaselist.in stationist.in y eventlist.in que son entrada para xcorloc. Los archivos se generan al hacer una consulta
a traves Seiscomp con parámetros ingresados por el usuario. 

Como entrada, la rutina requiere los rangos de tiempos, area de localizacion de los eventos de interés 
(p: poligono definido por vertices, c:circular o r:rectangular) magnitudes y errores. 
El ID de cada sismo se construye a partir del tiempo de origen (ver descripción de phase.dat en manual de hypoDD).
La rutina permite filtrar las estaciones de interés (que deben escribirse como una lista ["XXX","YYY",....] ), 
en este caso los archivos phase.dat  y station.dat omitirán las otras estaciones donde se pico el evento.


'''
import MySQLdb
import datetime
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.stream import Stream
import obspy
from obspy import read
import os
import shutil
from io import open
from matplotlib.dates import DateFormatter, date2num
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point, Polygon
import utm
from dateutil.relativedelta import relativedelta

#++++++-------------------------- parametros de usuario
time_i="2020-03-24 00:00:00"
time_f="2020-09-05 00:00:00"
#--------------parametros a filtrar
prof_min=-99
prof_max=50
mag_min=0
mag_max=10
elat_max=1000
elon_max=1000
eprof_max=2000
estaciones_interes=["CBOC","URE","ZAR", "MEDEC"] #si se desea generar el phase dat para ciertas estaciones, ponerlas como una lista y activar filtrar_estaciones con True
filtrar_estaciones=False


#-----------elegir tipo de busqueda
#se puede hacer tres tipos de busquedas espaciales de sismos para generar el phasedat c, r o p, en el caso de elegir p se debe cargar el poligno
#dependiendo de la busqueda seleccionar parametros rectangulares o circulares
tipo_busqueda="r" #elegir c: circular, r: rectangular, p: poligno

#busqueda rec
if tipo_busqueda=="r":
	lat_min=6.43
	lat_max=8.43
	lon_min=-74.7
	lon_max=-73.7

	#
	lat_centro=999 #pone estos valores por defecto porque la funcion los pide pero no los utiliza
	lon_centro=999
	radio=999 #en kilometros

#busqueda cir
if tipo_busqueda=="c":
	lat_centro=3.448
	lon_centro=-74.190
	radio=30 #en kilometros

	#
	lat_min=666 #pone estos valores por defecto porque la funcion los pide pero no los utiliza
	lat_max=999
	lon_min=-777
	lon_max=-666

#busqueda poligono en este caso ajusta automaticamente las coordendas de busqueda
#deben darse los vertices cada linea en un archivo txt (lon lat) la separacion entre lon y lat debe ser una tabulacion
if tipo_busqueda=="p":
	archivo_poligono="/home/mlizarazo/Documentos/codigos_python/seiscomp2phasedat/poligonos/poli_alg_norte.txt"
	arc=open(archivo_poligono,"r")
	lineas=arc.readlines()    
	coords=[]
	lons_poli=[]
	lats_poli=[]
	for l in lineas:
		
		vec_l=l.split("\t")
		coords.append((float(vec_l[0]), float(vec_l[1])))
		lons_poli.append(float(vec_l[0]))
		lats_poli.append(float(vec_l[1]))
	poligono = Polygon(coords)

	#define coordendas de busqueda en seicomp a traves de los maximos y minimos del poligno (a esa busqueda posteriormente se le filtrara en el polgno como tal)
	lat_min=min(lats_poli)
	lat_max=max(lats_poli)
	lon_min=min(lons_poli)
	lon_max=max(lons_poli)

	lat_centro=999 #pone estos valores por defecto porque la funcion los pide pero no los utiliza
	lon_centro=999
	radio=999 #en kilometros


#------Inicio de ejecucion
#-----------------------------------------------------------------

def escibe_phasedat(contador_sismos,c,a,b,fila,to,lat,lon,est_seiscomp,filas2,filtrar_estaciones,estaciones_interes,vec_amd,vec_hms,prof,mag,error_NE,error_Z,RMS,peso,archivo,archivo_xc,estaciones_loc,fecha,archivo_evenlist):

	#HYPODD soporta maximo 9 caracteres para el ID. 
	#Al intentar vincular la fecha como ID al ponerlo por ejmplo como yyyymmddhhmmss se tienen 14 caracteres 
	#asi mismo al convertirlo a dia juliano o date2num se obtienen numeros racionales con multiples decimales dada la precision de la fecha en segundos, 
	#con lo cual es imposible recuparear un numero de nueve digitos que represente la fecha
	#la solucion es generar un timestamp. Empricamente al restarle 30 años a la fecha actual daran timestamp con 8 y 9 digitos entre 
	#el año 2000 y 2031 (fechas validas para este codigo)
	#Como la idea es en algun caso recontruir la fecha de origen a partir del ID, para retornar a la fecha se debe ejecutar la siguiente linea
	# reconstruccion_fecha=datetime.datetime.fromtimestamp(int(ID_fecha)) + relativedelta(years=30)
	# En caso de requerir fechas anteriores al 2000 se generan timestamp negativos lo cual genera un problema para el ID, en este caso se deben restr mas años (mas de 30)
	# y se recomienda imprimir el ID_fecha hasta que no de negativo

	utc_fecha=UTCDateTime(fecha)
	t=utc_fecha.datetime+relativedelta(years=-30) #le quita a la fecha 1990 años, para que el ID pueda dsicriminarse y generarse con 7 digitos (el id de hipoDD soporta 9 digts y el de xcoorlog 7 digits)
	ID_fecha=str(int(datetime.datetime.timestamp(t)))
	
	#reconstruccion_fecha=datetime.datetime.fromtimestamp(int(ID_fecha)) + relativedelta(years=30)
	#print(reconstruccion_fecha)
	#formatos usa el input de xcorloc
	fmt_hdr = "%10s %4d %02d %02d %02d %02d %6.3f %9.5f %11.5f %7.3f %5.2f"
	fmt_ph = "%-2s %-5s %1s %7.3f %7.2f"
	

	estacion=fila[10]
	if estacion in est_seiscomp: #pgunta si la estacion que encontro en fases su informacion esta en la tabla de estaciones de sesicomp, de no ser no escribe esa fase en el phasedat pues no estara en el sattion.dat
							
		tipo_fase=fila[11]
		red=fila[18]
		time_us=int(fila[19])/1000000 #tiempo en microsengos se usara para phaselist
		sec_xcorloc = int(vec_hms[2])+time_us
		t_arr = UTCDateTime(fila[12]+datetime.timedelta(milliseconds=float(fila[13])/1000))
		t_absoluto_arr=round(t_arr-UTCDateTime(to),3)

		#-------- halla distancia epicentral.. input de phaselist.in
		x_evento = utm.from_latlon(lat, lon)[0]/1000#la coordenda en x es el elemento 0 de utm.from_latlon(lat, lon)
		y_evento = utm.from_latlon(lat, lon)[1]/1000 #en km

		index=est_seiscomp.index(estacion) #busca el indice en la lista est_estacion, los indices son los mismos que en filas2, salvo que en filas2 estan las coordendas
		#con ese incide entra a la fila de la estacion y obtiene [1] latitud y [2] longitud de la estacion
		lat_estacion=filas2[index][1]
		lon_estacion=filas2[index][2]
		elevacion_estacion=filas2[index][3]

		x_estacion=utm.from_latlon(lat_estacion,lon_estacion)[0]/1000
		y_estacion=utm.from_latlon(lat_estacion,lon_estacion)[1]/1000

		dist_epi=round( ( ((x_evento - x_estacion)**2)+((y_evento - y_estacion)**2) )**0.5, 2)
		#---------

		
		if t_absoluto_arr>0:

			if c==1:#--------------						
				
				if filtrar_estaciones==True:

					if estacion in estaciones_interes:
						contador_sismos=contador_sismos+1

						b=b+1

						archivo.writelines("#"+" "+vec_amd[0]+" "+vec_amd[1]+" "+vec_amd[2]+" "+vec_hms[0]+" "+vec_hms[1]+" "+vec_hms[2]+" "+str(lat)+" "+str(lon)+" "+str(prof)+" "+str(mag)+" "+str(error_NE)+" "+str(error_Z)+" "+str(RMS)+" "+ID_fecha+"\n")
						archivo.writelines(estacion+"\t"+str(t_absoluto_arr)+"\t"+str(peso)+"\t"+tipo_fase+"\n")
						
						archivo_xc.writelines(fmt_hdr %(ID_fecha,int(vec_amd[0]),int(vec_amd[1]),int(vec_amd[2]),int(vec_hms[0]),int(vec_hms[1]),sec_xcorloc,lat,lon,prof,mag)+"\n")
						archivo_xc.writelines(fmt_ph %(red,estacion,tipo_fase,t_absoluto_arr,dist_epi)+"\n")
												
						print("#"+" "+vec_amd[0]+" "+vec_amd[1]+" "+vec_amd[2]+" "+vec_hms[0]+" "+vec_hms[1]+" "+vec_hms[2]+" "+str(lat)+" "+str(lon)+" "+str(prof)+" "+str(mag)+" "+str(error_NE)+" "+str(error_Z)+" "+str(RMS)+" "+ID_fecha)
						archivo_evenlist.writelines(fmt_hdr %(ID_fecha,int(vec_amd[0]),int(vec_amd[1]),int(vec_amd[2]),int(vec_hms[0]),int(vec_hms[1]),sec_xcorloc,lat,lon,prof,mag)+"\n")
						
						estaciones_loc.append([estacion,red,lat_estacion,lon_estacion,elevacion_estacion])
				else:
					contador_sismos=contador_sismos+1
					archivo.writelines("#"+" "+vec_amd[0]+" "+vec_amd[1]+" "+vec_amd[2]+" "+vec_hms[0]+" "+vec_hms[1]+" "+vec_hms[2]+" "+str(lat)+" "+str(lon)+" "+str(prof)+" "+str(mag)+" "+str(error_NE)+" "+str(error_Z)+" "+str(RMS)+" "+ID_fecha+"\n")
					archivo.writelines(estacion+"\t"+str(t_absoluto_arr)+"\t"+str(peso)+"\t"+tipo_fase+"\n")
					
					archivo_xc.writelines(fmt_hdr %(ID_fecha,int(vec_amd[0]),int(vec_amd[1]),int(vec_amd[2]),int(vec_hms[0]),int(vec_hms[1]),sec_xcorloc,lat,lon,prof,mag)+"\n")
					archivo_xc.writelines(fmt_ph %(red,estacion,tipo_fase,t_absoluto_arr,dist_epi)+"\n")
						
					print("#"+" "+vec_amd[0]+" "+vec_amd[1]+" "+vec_amd[2]+" "+vec_hms[0]+" "+vec_hms[1]+" "+vec_hms[2]+" "+str(lat)+" "+str(lon)+" "+str(prof)+" "+str(mag)+" "+str(error_NE)+" "+str(error_Z)+" "+str(RMS)+" "+ID_fecha)
					archivo_evenlist.writelines(fmt_hdr %(ID_fecha,int(vec_amd[0]),int(vec_amd[1]),int(vec_amd[2]),int(vec_hms[0]),int(vec_hms[1]),sec_xcorloc,lat,lon,prof,mag)+"\n")
						

					estaciones_loc.append([estacion,red,lat_estacion,lon_estacion,elevacion_estacion])
					
			if c>1:#---------------------------
				if filtrar_estaciones==True:
					if estacion in estaciones_interes:

						#a y b sirven para rectificar si ya escribio la linea de localizacion cuando el usuario decide filtrar las estaciones, inicialmente c es un contador que encuntra en filas un ID particular
					#(filas tiene repetido el ID por cuantas fases haya de un sismo), cuando c==1 escribe la localizacion y cuando c>1 escribe las fases. Cuando pasa a otro ID en el for de IDs c vuelve a ser cero
					# !esto sirve para cuando no se filtran estaciones!. Pero cuando se filtran las estaciones, no necesariamente la primera vez que aparece el ID es una estacion de interes en la lista filas (generando problemas)
					#Por ello, se hizo el contador b y a
					#si b = 0 es porque no reconocio en la primera fila (del vector filas) (coicidente con el ID) una estacion de interes, asi que ese mensaje se envia hasta este if (en el que c>1), con el fin de reconocer 
					#que no se ha escrito la linea de localizacion, sin embargo, si solo se pone el condicional
					#de que b=0 para que escriba la linea de localizacion, para las demas estaciones que sean coincidentes con las estaciones de interes si lo escribiria, asi que el condicional a sirve para establecer cuantas 
					#veces encontro estaciones de interes (si encuentra mas de una estacion de interes no debe escribir la linea de localizacion ese numero de veces)
					#si a=1 es porque ya encontro al menos una estacion conicidente con las de interes, en ese caso escribe la linea de localizacion si y solo si b tambien es cero y si a=1, si no se cumple que a=1 y b=0 no escribe
					#linea de localizacion.
						a=a+1
						if b==0 and a==1:
							archivo.writelines("#"+" "+vec_amd[0]+" "+vec_amd[1]+" "+vec_amd[2]+" "+vec_hms[0]+" "+vec_hms[1]+" "+vec_hms[2]+" "+str(lat)+" "+str(lon)+" "+str(prof)+" "+str(mag)+" "+str(error_NE)+" "+str(error_Z)+" "+str(RMS)+" "+ID_fecha+"\n")
							
							archivo_xc.writelines(fmt_hdr %(ID_fecha,int(vec_amd[0]),int(vec_amd[1]),int(vec_amd[2]),int(vec_hms[0]),int(vec_hms[1]),sec_xcorloc,lat,lon,prof,mag)+"\n")
							
						archivo.writelines(estacion+"\t"+str(t_absoluto_arr)+"\t"+str(peso)+"\t"+tipo_fase+"\n")
						
						archivo_xc.writelines(fmt_ph %(red,estacion,tipo_fase,t_absoluto_arr,dist_epi)+"\n")
						
						estaciones_loc.append([estacion,red,lat_estacion,lon_estacion,elevacion_estacion])
				else:
					archivo.writelines(estacion+"\t"+str(t_absoluto_arr)+"\t"+str(peso)+"\t"+tipo_fase+"\n")
					
					archivo_xc.writelines(fmt_ph %(red,estacion,tipo_fase,t_absoluto_arr,dist_epi)+"\n")					
					
					estaciones_loc.append([estacion,red,lat_estacion,lon_estacion,elevacion_estacion])
				

		if t_absoluto_arr<0:
			print("\nRevisar este sismo, !tiene tiempos de arrivo negativos!")
			print("#"+" "+vec_amd[0]+" "+vec_amd[1]+" "+vec_amd[2]+" "+vec_hms[0]+" "+vec_hms[1]+" "+vec_hms[2]+" "+str(lat)+" "+str(lon)+" "+str(prof)+" "+str(mag)+" "+str(error_NE)+" "+str(error_Z)+" "+str(RMS)+" "+ID_fecha)
			print(estacion+"\t"+str(t_absoluto_arr)+"\t"+str(peso)+"\t"+tipo_fase+"\n")


	return contador_sismos,c,b,a,archivo, archivo_xc, estaciones_loc,archivo_evenlist

def depurador_de_fases(filas,filas2):
	'''
	#----------------corrige error de fases repetidas en eventos anteriores a 01/03/2019 --------------------

	Depura algunas fases que estan erradas para los eventos que estuvieron en el proceso de la migracion, en el que hay fases repetidas o no se distingue P de S, 
	si se activa la depuracion imprimira cual fue el cambio, generalmente no encuentra problema en las fases
	La depuracion consiste en que si el t origen, la estacion, la magnitud y la fase son las mismas entra a preguntar
	si el tiempo absoluto tambien es el mismo, de serlo se trata de un evento repetido, pero si no lo es, es porque en la migracion se cambio 
	la fase (P a S o visceverza) en este caso para el tiempo absoluto menor adjudica la fase P y para el mayor la fase S
	En la consulta se cambio Pick.phaseHint_code por Arrival.phase_code, para no repitir las fases, sin embargo, tras la migracion algunos 
	fases aparecen doblemente picadas sin que se deba a un probelma humano. De esta manera el siguiente apartado de codigo calcula un ts - tp teorico 
	en una sola capa, con una velocidad de S de 3.5 km/s y una relacion Vp = raiz(3)Vs.
	Dicho S-P teorico es infraestimado (concientemente) dada la suposicion de una sola capa. De esta manera si se encuentra que la diferencia del tiempo 
	de arribo entre dos fases es menor a ese delta, no corresponde a un problema de mal etiquetada cierta fase (P por S) sino mas bien una picada doble 
	de la misma fase (error que puede provenir desde el sfile de Seisan) con ello debe ser eliminada una de esas fases. De esta forma 
	si hay dos fases repetidas en la misma estacion, con el mismo tiempo de origen, la misma localizacion y su delta observado es menor al delta teorico 
	entran a revision para eliminar una de ellas; de lo contrayio si una de las fases tiene polaridad y la otra no se asume que la primera es la correcta 
	y se elimina la otra, si no tienen polaridad ninguna de las dos se elimina la de mayor residual entre ambas. Si todo es igual incluyendo
	residual se elimina una de ellas al azar.
	'''
	
	e=0
	while e < len(filas)-1:

		e=e+1
		filas[e-1]=list(filas[e-1])
		filas[e]=list(filas[e])			
		
		if filas[e][1] == filas[e-1][1] and filas[e][10]==filas[e-1][10] and filas[e][5] == filas[e-1][5] and filas[e][11] == filas[e-1][11]: 
						
			for est in filas2: #filas2 se refiere al resultado de la conulta de todas las estaciones en seiscomp
				nombre_est=est[0]
				if nombre_est in filas[e][10]:
					lat_est=round(float(est[1]),4)
					lon_est=round(float(est[2]),4)
			'''
			Hasta aqui encontro el mismo evento picado en la misma estacion donde esta repetida la etiqueta de la fase P o la S.
			para determinar de cual se trata (antes de eliminarlas) calcula en una sola capa de baja velocidad el tiempo teorico S - P 
			si S-P teorico es mayor que la diferencia de tiempos abosultos observados de las dos fases en cuestion
			entra a decidir cual eliminar basado en residuales, polaridad de ambas fases. Por el contratio Si se cumple que el tiempo 
			teorico da para interpretar que se trata de S y P, hace el cambio de la etiqueta de la fase

			'''
			r = ( ((abs(filas[e][2]-lat_est)*111)**2) + ((abs(filas[e][3]-lon_est)*111)**2) + (filas[e][4]**2) )**0.5 
			del_t_teorico= 0.121*r

			if abs(UTCDateTime(filas[e][12]+datetime.timedelta(milliseconds=float(filas[e][13])/1000)) - UTCDateTime(filas[e-1][12]+datetime.timedelta(milliseconds=float(filas[e-1][13])/1000))) < del_t_teorico: 
						
				print("\n--------fases repetidas\n")					
				print(filas[e-1])
				print(filas[e])
				print("distancia teorica (una sola capa): "+str(r)+"km")
				print("ts  - tp teorico: "+str(del_t_teorico)+"s")	
				print("\ndelta de tiempo "+str(abs(UTCDateTime(filas[e][12]+datetime.timedelta(milliseconds=float(filas[e][13])/1000)) - UTCDateTime(filas[e-1][12]+datetime.timedelta(milliseconds=float(filas[e-1][13])/1000)))))

				polaridad1=filas[e-1][16]
				polaridad2=filas[e][16]
				residual1=filas[e-1][17]
				residual2=filas[e][17]

				#-------------criterios de eliminacion

				if polaridad1==None and polaridad2!=None:					

					print("\nse elimino pues no tiene polaridad: \n"+str(filas[e-1])+"\n")
					filas.pop(e-1) #con el indice elimina un elemento de una lista
				
				if polaridad1!=None and polaridad2==None:					

					print("\nse elimino pues no tiene polaridad: \n"+str(filas[e])+"\n")
					filas.pop(e) #con el indice elimina un elemento de una lista

				if polaridad1==None and polaridad2==None:					

					if residual1==None and residual2!=None:
						print("\nse elimino pues no tiene residual ni polaridad: \n"+str(filas[e-1])+"\n")
						filas.pop(e-1) #con el indice elimina un elemento de una lista

					if residual1!=None and residual2==None:
						print("\nse elimino pues no tiene residual ni polaridad: \n"+str(filas[e])+"\n")
						filas.pop(e) #con el indice elimina un elemento de una lista

					if residual1==None and residual2==None:
						print("\nse elimino cualquiera de los dos pues no tienen residual ni polaridad: \n"+str(filas[e])+"\n")
						filas.pop(e) #con el indice elimina un elemento de una lista
						#se elimina cualquiera

					if residual1!=None and residual2!=None:
						
						residual1=abs(residual1)
						residual2=abs(residual2)

						if residual1 > residual2:

							print("\nse elimino pues tiene mayor residual: \n"+str(filas[e-1])+"\n")
							filas.pop(e-1) #con el indice elimina un elemento de una lista
							

						if residual1 < residual2:

							print("\nse elimino pues tiene mayor residual: \n"+str(filas[e])+"\n")
							filas.pop(e) #con el indice elimina un elemento de una lista
							

						if residual1 == residual2:

							print("\nse elimino cualquiera de los dos pues no tienen polaridad y tienen el mismo residual: \n"+str(filas[e])+"\n")
							filas.pop(e) #con el indice elimina un elemento de una lista
							#se elimina cualquiera

				if polaridad1!=None and polaridad2!=None:					

					if residual1==None and residual2!=None:
						print("\nse elimino pues no tiene residual: \n"+str(filas[e-1])+"\n")
						filas.pop(e-1) #con el indice elimina un elemento de una lista

					if residual1!=None and residual2==None:
						print("\nse elimino pues no tiene residual: \n"+str(filas[e])+"\n")
						filas.pop(e) #con el indice elimina un elemento de una lista

					if residual1==None and residual2==None:
						print("\nse elimino cualquiera de los dos pues no tienen residual: \n"+str(filas[e])+"\n")
						filas.pop(e) #con el indice elimina un elemento de una lista
						#se elimina cualquiera

					if residual1!=None and residual2!=None:
						
						residual1=abs(residual1)
						residual2=abs(residual2)

						if residual1 > residual2:

							print("\nse elimino pues tiene mayor residual: \n"+str(filas[e-1])+"\n")
							filas.pop(e-1) #con el indice elimina un elemento de una lista
							

						if residual1 < residual2:

							print("\nse elimino pues tiene mayor residual: \n"+str(filas[e])+"\n")
							filas.pop(e) #con el indice elimina un elemento de una lista
							

						if residual1 == residual2:

							print("\nse elimino cualquiera de los dos pues tienen el mismo residual: \n"+str(filas[e])+"\n")
							filas.pop(e) #con el indice elimina un elemento de una lista
							#se elimina cualquiera
				#---------------------

			#en este else habria reconocido una fase mal etiquetada y se reemplaza P por S o visceverza
			else:
				print("este evento tiene las fases mal etiquetadas o repetidas (por ejemplo en la base de datos aparecen dos S picadas en N y E), por defecto se cambiara S por P o visceverza, sin embargo revisar!! ")
				print(filas[e])
				print(filas[e-1])

				if UTCDateTime(filas[e][12]+datetime.timedelta(milliseconds=float(filas[e][13])/1000)) > UTCDateTime(filas[e-1][12]+datetime.timedelta(milliseconds=float(filas[e-1][13])/1000)):
					filas[e][11] = "S"
					filas[e-1][11] = "P"
				else:
					filas[e][11] = "P"
					filas[e-1][11] = "S"
				
				print("se reemplazaron las fases del evento "+ str(filas[e][1]))
				print(filas[e])
				print(filas[e-1])

	return filas

def seiscomp2phasedat(time_i,time_f,lat_min,lat_max,lon_min,lon_max,lat_centro,lon_centro,radio,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,estaciones_interes,filtrar_estaciones,tipo_busqueda):


	if tipo_busqueda=="c":
		
		epsil=(2*radio/111)
		mapa = Basemap(lon_centro-epsil,lat_centro-epsil,lon_centro+epsil,lat_centro+epsil, epsg=3117, resolution=None) #grafica el rectangulo de lat y lon max y min
		x_centro,y_centro=mapa(lon_centro,lat_centro)
		
		#modifica los minimos que ingresan a la consulta que por defecto es rectangular
		lat_min=lat_centro-epsil
		lat_max=lat_centro+epsil
		lon_min=lon_centro-epsil
		lon_max=lon_centro+epsil


	#db = MySQLdb.connect(host="10.100.100.232",user="consulta", passwd="consulta", db="seiscomp3")
	db = MySQLdb.connect(host="172.25.3.135", user="consulta", passwd="consulta", db="seiscomp3")   

	cur = db.cursor()

	cur.execute("Select POEv.publicID, \
	Origin.time_value as 'Fecha', \
	ROUND(Origin.latitude_value,6) as 'Lat', \
	ROUND(Origin.longitude_value,6) as 'Lon', \
	ROUND(Origin.depth_value,2) as 'Z', \
	ROUND(Magnitude.magnitude_value,1) as 'Mag', \
	ROUND(Origin.latitude_uncertainty,3), \
	ROUND(Origin.longitude_uncertainty,3), \
	ROUND(Origin.depth_uncertainty,2), \
	ROUND(Origin.quality_standardError,1) as 'RMS', \
	Pick.waveformID_stationCode, \
	Arrival.phase_code, \
	Pick.time_value, \
	Pick.time_value_ms, \
	Arrival.weight, \
	Pick.evaluationMode, \
	Pick.polarity, \
	Arrival.timeResidual, \
	Pick.waveformID_networkCode, \
	Origin.time_value_ms \
	from Event AS EvMF \
	left join PublicObject AS POEv ON EvMF._oid = POEv._oid \
	left join PublicObject as POOri ON EvMF.preferredOriginID=POOri.publicID \
	left join Origin ON POOri._oid=Origin._oid \
	left join PublicObject as POMag on EvMF.preferredMagnitudeID=POMag.publicID \
	left join Magnitude ON Magnitude._oid = POMag._oid \
	left join Arrival on Arrival._parent_oid=Origin._oid \
	left join PublicObject as POOri1 on POOri1.publicID = Arrival.pickID \
	left join Pick on Pick._oid= POOri1._oid \
	left join EventDescription ON EvMF._oid=EventDescription._parent_oid \
	where \
	(Magnitude.magnitude_value between "+str(mag_min)+" and "+str(mag_max)+" and \
	Origin.latitude_value between "+str(lat_min)+" and "+str(lat_max)+" and \
	Origin.longitude_value between "+str(lon_min)+" and "+str(lon_max)+" and \
	Origin.depth_value between "+str(prof_min)+" and "+str(prof_max)+" and \
	Origin.latitude_uncertainty<"+str(elat_max)+" and \
	Origin.longitude_uncertainty<"+str(elon_max)+" and \
	Origin.depth_uncertainty<"+str(eprof_max)+" and \
	Origin.time_value between '"+time_i+"' and '"+time_f+"') and \
	Pick.evaluationMode = 'manual' and \
	(EvMF.type like 'earthquake' or EvMF.type is null )") #ultima linea para filtrar los no localizables (para casos como el de providencia quitar esta linea, pues interesan los no localizables)

	filas = cur.fetchall()
	filas=list(filas)

	#--------------
	cur2 = db.cursor()
	
	#hace una consulta de todas las estaciones de la red, lo cruza con las estaciones almacenadas arriba en estaciones_loc y con ello escribe las coordenadas para cada una de esas estaciones en el archivo station.dat
	cur2.execute("Select Station.code, \
				Station.latitude, \
				Station.longitude, \
				Station.elevation \
				from Station") 
	
	filas2 = cur2.fetchall()
	filas2=list(filas2) #filas2 es la lista de todas las estaciones en seiscomp con sus coordenadas: name_est, lat, lon

	est_seiscomp=[]
	for est in filas2: #filas2 se refiere al resultado de la conulta de todas las estaciones en seiscomp con sus coordenadas: name_est, lat, lon
		est_seiscomp.append(est[0])
	#--------------------			

	archivo=open("phase.dat","w") #escribira el archivo phase.dat input de ph2dt/hypoDD
	archivo_xc=open("phaselist.in","w") #escribira el archivo phaselist.in input de xcoorlog
	archivo_evenlist=open("eventlist.in","w") #escribira listado de eventos input de xcoorlog

	#----------------corrige error de fases repetidas en eventos anteriores a 01/03/2019 --------------------
	filas=depurador_de_fases(filas,filas2)
	#-------------------------
	
	#------------------------genera vector con IDs unicos y ordenados
	IDs=[]
	for i in range(1,len(filas)):
		ID1=filas[i-1][0]
		ID2=filas[i][0]
		if ID1!=ID2:
			IDs.append(ID1)
	IDs.append(filas[-1][0]) 
	#--------------------------------------------------------------
	print("\nSismos encontrados: \n")

	#-------------------escribe archivo phase.dat y sattion.dat
	estaciones_loc=[] #almacenara las estaciones en las que se registraron los sismos encontrados
	
	contador_sismos=0
	for ID in IDs:
		
		c=0
		b=0
		a=0

		for fila in filas:
			
			if fila[0]==ID:
				
				c=c+1
				to=fila[1]
				#print(to)
				fecha=str(to)
				vec_fecha=fecha.split(" ")
				an_me_di=vec_fecha[0]
				ho_mi_se=vec_fecha[1]
				vec_amd=an_me_di.split("-")
				vec_hms=ho_mi_se.split(":")
				lat=fila[2]
				lon=fila[3]
				prof=fila[4]
				mag=fila[5]
				error_NE=fila[6]
				error_Z=fila[8]
				RMS=fila[9]
				peso=round(float(fila[14]),1)
				#----- correccion del peso: ojo en la migracion vienen pesos de valores hasta 50, 100 quizas porque el localizador no era hypo71 como en seiscomp (que inicio en 2018)
				#como se migraron desde 2016 quedaron revueltos los pesos...  en hypoDD el maximo es 1 entonces hay que dividir en 100 a aquellos de esa escala 
				if peso == 100.0 or peso == 75.0 or peso == 50.0 or peso == 25.0:
					peso=round(peso/100,1)
					

				if tipo_busqueda=="c":
			
					x_sis,y_sis=mapa(lon,lat)			

					if (x_centro - (radio*1000)) < x_sis < (x_centro + (radio*1000)):
				
						if  y_centro - (((radio*1000)**2)-(x_sis-x_centro)**2)**(0.5) < y_sis < y_centro + (((radio*1000)**2)-(x_sis-x_centro)**2)**(0.5):

							contador_sismos,c,b,a,archivo_xc,archivo_xc,estaciones_loc,archivo_evenlist=escibe_phasedat(contador_sismos,c,a,b,fila,to,lat,lon,est_seiscomp,filas2,filtrar_estaciones,estaciones_interes,vec_amd,vec_hms,prof,mag,error_NE,error_Z,RMS,peso,archivo,archivo_xc,estaciones_loc,fecha,archivo_evenlist)
				
				if tipo_busqueda=="r":

					contador_sismos,c,b,a,archivo_xc,archivo_xc,estaciones_loc,archivo_evenlist=escibe_phasedat(contador_sismos,c,a,b,fila,to,lat,lon,est_seiscomp,filas2,filtrar_estaciones,estaciones_interes,vec_amd,vec_hms,prof,mag,error_NE,error_Z,RMS,peso,archivo,archivo_xc,estaciones_loc,fecha,archivo_evenlist)

				if tipo_busqueda=="p":

					p = Point(lon, lat)

					if p.within(poligono)==True:

						contador_sismos,c,b,a,archivo_xc,archivo_xc,estaciones_loc,archivo_evenlist=escibe_phasedat(contador_sismos,c,a,b,fila,to,lat,lon,est_seiscomp,filas2,filtrar_estaciones,estaciones_interes,vec_amd,vec_hms,prof,mag,error_NE,error_Z,RMS,peso,archivo,archivo_xc,estaciones_loc,fecha,archivo_evenlist)


	archivo.close()
	archivo_xc.close()
	archivo_evenlist.close()

	print("Se encontraron "+ str(contador_sismos)+" sismos")

	
	#-------------------------escribe el station.dat 
	archivo_stationdat=open("station.dat","w") #escribira el station.dat para ph2dt/hypodd
	archivo_stationlist=open("stationlist.in","w") #escribira el stationlist.in para Xcorloc

	set_est=[]	
	#estaciones_loc contiene [estacion,red,lat_estacion,lon_estacion,elevacion_estacion] (estan repetidas las estaciones por cuantas fases tiene un evento y numero de eventos)
	
	for est in estaciones_loc:
		nombre_est=est[0]
		if nombre_est not in set_est:
			red=est[1]
			lat_est=round(float(est[2]),4)
			lon_est=round(float(est[3]),4)
			elev_est=round(float(est[4]),3)
			archivo_stationdat.writelines(nombre_est+"\t"+str(lat_est)+"\t"+str(lon_est)+"\n")
			archivo_stationlist.writelines(red+"\t"+nombre_est+"\t"+str(lat_est)+"\t"+str(lon_est)+"\t"+str(elev_est)+"\n")
		set_est.append(nombre_est)


	archivo_stationdat.close()
	archivo_stationlist.close()
	print("\nListo, !se guardaron en el directorio los archivos phase.dat, station.dat (input de ph2dt/hypoDD), phaselist.in y stlist.in (input de xcorloc)")

seiscomp2phasedat(time_i,time_f,lat_min,lat_max,lon_min,lon_max,lat_centro,lon_centro,radio,prof_min,prof_max,mag_min,mag_max,elat_max,elon_max,eprof_max,estaciones_interes,filtrar_estaciones,tipo_busqueda)

