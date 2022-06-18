/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: ParedMCPS
 CREATION DATE: 20/04/2022
-----------------------------------------------------------------------------------------*/
USE MATH
COMPONENT Pared
	PORTS
		OUT SRB salida
	DATA
		REAL Tc = 3330
		REAL D0 = 0.0254 UNITS u_m   "Diametro inicial geometria cilindrica"
	INIT			----------Configuracion geometrica inicial, esto es provisional, se tendria que introducir en el experimento?-----------
			salida.Y = 0
			salida.S = 2*MATH.PI*0.5*D0
		   salida.Ap = MATH.PI*(0.5*D0)**2
	CONTINUOUS 
	salida.Tt = Tc		
	salida.Ucin = 0
	salida.Apcin = salida.Ap
	salida.dgin = 0
	salida.g = 0
	salida.gprev = 0
	salida.Coord = 0
	salida.Combustion = 1
	
	
	salida.Pout = salida.Pin
	salida.r0in= salida.r0out
	
	
END COMPONENT