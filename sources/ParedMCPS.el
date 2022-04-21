/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: ParedMCPS
 CREATION DATE: 20/04/2022
-----------------------------------------------------------------------------------------*/
COMPONENT Pared
	PORTS
		OUT test salida
	DATA
		REAL P = 5E6
		REAL U = 0 
		REAL Coord = 0  UNITS "m"	
		REAL g = 0

	DECLS
	
	CONTINUOUS 
	
	--P = 5E6
	salida.Pin = P
	salida.Pout = P
	salida.U = U
	salida.Coord = Coord
	salida.g = g
	
	
END COMPONENT