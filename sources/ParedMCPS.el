/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: ParedMCPS
 CREATION DATE: 20/04/2022
-----------------------------------------------------------------------------------------*/
COMPONENT Pared
	PORTS
		OUT test salida
	DATA
		REAL Tc = 3330
		
	DECLS
		BOUND REAL Pheadend
	INIT
	
	CONTINUOUS 
	
	salida.Pout = salida.Pin
	salida.Pin = Pheadend
	salida.Tin = salida.Tout
	salida.Tout = Tc		
	salida.U = 0
	salida.g = 0
	salida.Coord = 0
	
	
	
END COMPONENT