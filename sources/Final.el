/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: Final
 CREATION DATE: 20/04/2022
-----------------------------------------------------------------------------------------*/
USE LPRES

COMPONENT Conversor
	PORTS
		IN test entrada
		--OUT Fluid salida
	DATA
		REAL c_star = 1500
		REAL A_th = 8.00E-03
	DECLS
		BOUND REAL g 	
		--BOUND REAL Ptdischarge
		--BOUND REAL Pdischarge
	CONTINUOUS		
		entrada.Pin = entrada.Pout
		--Pdischarge = entrada.Pout
		
		
	
		-- Condicion de igualdad de gastos entre combustion y tobera 
		--1 = (g*c_star)/(A_th*Ptdischarge)  --Habria que imponerla en los boundary (presion de remanso en el inicio)-- 
		
		c_star=entrada.Pin*A_th/entrada.g

		-- Hay que especificar el gas de salida para que funcione la tobera
		--ENUM Gases salida.Fluid = setofElem(gases,1)		
		

END COMPONENT