/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: ParedMCPS
 CREATION DATE: 20/04/2022
-----------------------------------------------------------------------------------------*/
COMPONENT Pared
	PORTS
		OUT test salida
	DATA
	
	DECLS
	
	CONTINUOUS 
	
	salida.Pout = salida.Pin
	salida.Tin = salida.Tout
	salida.Tout = 3300
	salida.U = 0
	salida.Coord = 0
	salida.g = 0
	
	
END COMPONENT