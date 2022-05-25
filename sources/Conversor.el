/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: Final
 CREATION DATE: 20/04/2022
-----------------------------------------------------------------------------------------*/
USE LPRES

COMPONENT Conversor
	PORTS
		IN test entrada
	DATA
		REAL c_star = 1500
		REAL A_th = 8.00E-03
		REAL R_gas = 307.93
		REAL gamma = 1.33
	DECLS
		BOUND REAL Ptdischarge
		REAL Rho
		REAL SoundSpeed
		REAL MACH
		--REAL T
		
	CONTINUOUS		
		
		entrada.Pin = entrada.Pout
		entrada.Tin = entrada.Tout
		Rho = entrada.Pin /(R_gas*entrada.Tin)
		SoundSpeed = sqrt(gamma*R_gas*entrada.Tin)
		MACH = entrada.U/SoundSpeed
		Ptdischarge = entrada.Pin  *(1+0.5*(gamma-1)*MACH**2)**(gamma/(gamma-1))
	
		-- Condicion de igualdad de gastos entre combustion y tobera 
		1 = ((entrada.g*c_star)/(A_th*Ptdischarge))  --Habria que imponerla en los boundary (presion de remanso en el inicio)-- 
	
		--T = entrada.g*Usalida

END COMPONENT