/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: PuertoMCHsolido
 CREATION DATE: 30/04/2022
-----------------------------------------------------------------------------------------*/
USE MATH VERSION "3.2"

PORT test SINGLE
		
SUM	                   REAL g                 UNITS u_kg_s         RANGE ZERO,Inf       "Mass flow"
EQUAL		                REAL Ap                UNITS u_m2           RANGE ZERO,Inf       "Port Area"
EQUAL							 REAL S					   UNITS u_m				RANGE ZERO,Inf			"Perimeter of port"
EQUAL		                REAL U                 UNITS u_m2           RANGE ZERO,Inf       "Flow speed"			
EQUAL							 REAL Coord				   UNITS u_m				RANGE ZERO,Inf			"Coordinate"

EQUAL							 REAL Pout				   UNITS u_m				RANGE ZERO,Inf		   "Pre-interface pressure"
EQUAL							 REAL Pin				   UNITS u_m				RANGE ZERO,Inf		   "Post-interface pressure"
EQUAL							 REAL Tout				   UNITS u_m				RANGE ZERO,Inf		   "Pre-interface temp"
EQUAL							 REAL Tin				   UNITS u_m				RANGE ZERO,Inf		   "Post-interface temp"	
-- Prueba
REAL rp0_intermedia						

CONTINUOUS
rp0_intermedia = (1.27955E-05)* (0.5*(Pin+Pout))**0.4    -- ¿Hay que pasar variables de la ley de vieille?
S'  =	2*MATH.PI*rp0_intermedia
Ap' = S*rp0_intermedia
END PORT