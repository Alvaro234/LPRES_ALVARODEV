/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: Final
 CREATION DATE: 20/04/2022
-----------------------------------------------------------------------------------------*/
USE LPRES
USE MATH
COMPONENT Conversor
	PORTS
		IN SRB entrada
	DATA
		-- Datos gas y propelente --
		REAL c_star = 1500
		REAL R_gas = 307.93
		REAL gamma = 1.33
		--- Geometria tobera ---
		REAL D = 0.0635 UNITS u_m	  "Diametro carcasa cilindrica"
		REAL D0 = 0.0254 UNITS u_m   "Diametro inicial geometria cilindrica"
		REAL A_th = 8.00E-03
	DECLS
		REAL Pt
		REAL Rho
		REAL SoundSpeed
		REAL MACH
	
	INIT			----------Configuracion geometrica inicial, esto es provisional, se tendria que introducir en el experimento?-----------	
		-- Habria que implementar una distribucion de presiones inicial a partir de P0 = (rho_p*a*c_estrella* Abtotal/At)**(1/(1-n))
		entrada.Y = 0
		entrada.S = 2*MATH.PI*0.5*D0
		entrada.Ap = MATH.PI*(0.5*D0)**2
	
	CONTINUOUS		
		

		entrada.Apcout = MATH.PI*(0.5*D0)**2    // Un area caracteristica de la "pretobera"
		entrada.Ucout = entrada.g/Rho*entrada.Apcout--entrada.Ucin     //Algo con la tobera
		--entrada.Pout = entrada.Pin  	//Perdida de presion de un 10% ¿?
		entrada.r0in= entrada.r0out
		

		Rho = entrada.Pout /(R_gas*entrada.Tt)
		SoundSpeed = sqrt(gamma*R_gas*entrada.Tt)
		MACH = entrada.Ucout/SoundSpeed
		Pt = entrada.Pout  *(1+0.5*(gamma-1)*MACH**2)**(gamma/(gamma-1))--Aqui se le podria poner un coeficiente de 0,9 de perdida de presion entre pastilla y tober
		
		
		
		------------------------- ECUACION DE CIERRE -------------------------------------------
		c_star = Pt*A_th/entrada.g   -- Condicion de igualdad de gastos entre combustion y tobera 
		
		
		
	
	
		 ------ Erosion cte de la garganta -----
		 /*
		 rth' = alfa			-- Valor inicial rth = D0/2
		 At=MATH.PI*(rth)**2
		 
		
		 T =  Ct*A_th*entrada.Pout
		 Ct = Fgamma*sqrt( (2*gamma/(gamma-1))*(1-(Pe/entrada.Pout)**((gamma-1)/gamma))  )
		 ER = 
		  */
		------------- ECUACIONES TOBERA  -----------------

		

END COMPONENT