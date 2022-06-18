 /*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: PuertoMCHsolido
 CREATION DATE: 30/04/2022
-----------------------------------------------------------------------------------------*/
USE MATH VERSION "3.2"

PORT SRB SINGLE
	
EQUAL		                REAL Ap                UNITS u_m2           RANGE ZERO,Inf       "Port Area"
EQUAL							 REAL S					   UNITS u_m				RANGE ZERO,Inf			"Perimeter of port"
EQUAL							 REAL Y					   UNITS u_m				RANGE ZERO,Inf			"Lateral coordinate"			
			
EQUAL							 REAL Coord				   UNITS u_m				RANGE ZERO,Inf			"Coordinate"
EQUAL							 REAL Tt				     UNITS u_K				RANGE ZERO,Inf		   "Total temperature"
EQUAL                    REAL g                 UNITS u_kg_s         RANGE ZERO,Inf       "Mass flow"	
--EQUAL		                REAL U                 UNITS u_m_s           RANGE ZERO,Inf       "Flow speed"			
						 		 
EQUAL							 REAL Pout				   UNITS u_Pa				RANGE ZERO,Inf		   "Pre-interface pressure" --mirar
EQUAL							 REAL Pin				   UNITS u_Pa				RANGE ZERO,Inf		   "Post-interface pressure"
EQUAL							 REAL dgout				   UNITS u_kg_s			RANGE ZERO,Inf		   "Pre-interface pressure"
EQUAL							 REAL dgin				   UNITS u_kg_s			RANGE ZERO,Inf		   "Post-interface pressure"
EQUAL							 REAL Ucout				   UNITS u_m_s				RANGE ZERO,Inf		   "Pre-interface pressure"
EQUAL							 REAL Ucin				   UNITS u_m_s				RANGE ZERO,Inf		   "Post-interface pressure"
EQUAL							 REAL Apcout				UNITS u_m2				RANGE ZERO,Inf		   "Pre-interface pressure"
EQUAL							 REAL Apcin				   UNITS u_m2				RANGE ZERO,Inf		   "Post-interface pressure"
EQUAL                    REAL gprev             UNITS u_kg_s         RANGE ZERO,Inf       "Mass flow"	

EQUAL							 REAL r0in				   UNITS u_m				RANGE ZERO,Inf		   "Pre-interface erosion"
EQUAL						    REAL r0out				   UNITS u_m				RANGE ZERO,Inf		   "Post-interface erosion"
EQUAL					    	 REAL Combustion				   UNITS u_m				RANGE ZERO,Inf "Combustion state"
		

CONTINUOUS
			Y' = 	0.5*(r0out+r0in)*Combustion  -- Revisar
			S'  =	2*MATH.PI*Y'     				
			Ap' = S * Y'

		  (g+0.5*dgout)*Ucout + Apcout*Pout = (gprev+0.5*dgin)*Ucin+Apcin*Pin  -- Ec cant mov

END PORT