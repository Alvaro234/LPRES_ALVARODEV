/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: Final
 CREATION DATE: 20/04/2022
-----------------------------------------------------------------------------------------*/
USE LPRES

COMPONENT Final
	PORTS
		IN test entrada
		--OUT Fluid salida
	DATA
		REAL c_star = 1500
		REAL A_th = 8.00E-03
		REAL pt = 5.83E+06

	DECLS
	
		REAL g 
		
					
	CONTINUOUS		
		-- Hay que especificar el gas de salida para que funcione la tobera
		--ENUM Gases salida.Fluid = setofElem(gases,1)		
		
		
		-- Condicion de igualdad de gastos entre combustion y tobera --  ¿Tramos intermedios? 
		--1 = (g*c_star)/(A_th*pt)  --Habria que imponerla en los boundary (presion de remanso en el inicio)-- 

END COMPONENT