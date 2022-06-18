/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: MChSolido2
 AUTHOR: Álvaro González Villarreal
 DESCRIptION: Seccion de motor cohete de propelente solido V0.2
 CREATION DATE: 02/04/2022
-----------------------------------------------------------------------------------------*/

USE MATH VERSION "3.2"
USE LPRES

COMPONENT MChSolido2(INTEGER N = 3, BOOLEAN Comb_Erosiva = TRUE )
	PORTS
		IN  SRB entrada
		OUT SRB salida
	DATA
			
			-- DATOS GENERALES COMBUSTIBLE Y PARAMETROS COMBUSTION
			REAL gamma = 1.33
			REAL Rho_P = 1750							 "Densidad del propelente"
			REAL M_m = 26.144    UNITS u_g 	 "Masa molar del propelente"
			REAL R_gas = 307.93						 "Constante R del gas"
			REAL Tc = 3300 							 "Temperatura de combustion adiabatica"
			-- Ley de vieille  
			REAL Exponente =	0.4					 "Exponente de Vieille"
			REAL Factor_a =	1.27955E-05			 "Factor de Vieille"	
			-- Combustion erosiva
			REAL Viscosidad_camara =	1.80E-05	 "Viscosidad de camara"
		   REAL gth =	35								 "Parametro de erosion"
			-- DATOS GEOMETRIA
			REAL L = 5 UNITS u_m 		  "Longitud de la seccion"
			REAL D = 0.0635 UNITS u_m	  "Diametro carcasa cilindrica"
			REAL D0 = 0.0254 UNITS u_m   "Diametro inicial geometria cilindrica"
	DECLS
				------ Interfaces (1 a N+1)----------------
				-- Parametros geometricos
				DISCR REAL dx	
				REAL Y[N+1]
				REAL S[N+1]
				REAL Ap[N+1]
				
				REAL Ab_total
				REAL V_total
				REAL Masa_total
				REAL LoadFraction
				
				
				
				REAL Coord [N+1]  UNITS "m"	
				REAL g[N+1]
				DISCR REAL Combustion[N+1]
				DISCR REAL A_carcasa
				---- "Volumenes" discretos	------ (1 a N)		
				REAL Ab[N] 
				REAL dg[N]
				REAL P[N]        UNITS u_Pa				RANGE ZERO,Inf		   "Pressure"      
				REAL T[N]			UNITS u_Pa				RANGE ZERO,Inf		   "Temperature"      
				REAL Uc[N]
				REAL Apc[N]
		      -----------------------
				REAL rp0[N]
				REAL r0[N]
				REAL Rho[N] 
				REAL SoundSpeed[N]
				REAL MACH[N]
				REAL Tt[N+1]
				REAL pt[N]
				
				
				-- Parametros de combustion erosiva --
				-- Modelo de mukunda --
				REAL g_i[N]
				REAL g0_i[N]  
				REAL Re0[N]
				REAL eta_temp[N]
				REAL eta[N]				RANGE 1,Inf
				
				
				
	INIT		-----------Configuracion geometrica inicial, esto es provisional, se tendria que introducir en el experimento-----------
		dx = L/N	-- Aqui se podria poner una pos diferente para cada interfaz y no una malla equiespaciada
		FOR(i IN 1,N+1)
			Y[i] = 0
			S[i] = 2*MATH.PI*0.5*D0
		   Ap[i] = MATH.PI*(0.5*D0)**2
			A_carcasa = MATH.PI*(D/2)**2
			Combustion[i] = 1
		END FOR
	DISCRETE
		EXPAND(i IN 1,N+1)
			WHEN (Y[i] >= (D-D0)/2) TOL 1E-05 THEN
				Combustion[i] = 0			// La combustion en una interfaz acaba cuando se llega a la carcasa
			END WHEN
			WHEN(SUM(i IN 1,N+1; Combustion[i]) <= (N+1)*0.6 ) THEN
		      --STOP "******************NO FUEL LEFT*********************" // Realmente es que queden 3 interfaces o menos por quemar
				// Se podria hacer cuando la masa de propelente sea un 1-5% de la inicial (reg stacionario)? ¿Como lo hacemos cuando haya varios componentes?
			END WHEN
	CONTINUOUS		
		----------------------------------------------------------------INTERFACES (De 1 a N+1)---------------------------------------------------------------------
		EXPAND_BLOCK(i IN 1,N)
			Coord[i+1] = Coord[i] + dx			 
			g[i+1] = g[i]+dg[i]
		END EXPAND_BLOCK	
		EXPAND_BLOCK(i IN 2,N)
			Y[i]' = 	0.5*(r0[i]+r0[i-1])*Combustion[i]  -- Revisar (parece funcionar) 
			S[i]'  =	2*MATH.PI*Y[i]'     -- La idea es poner una tabla para 1 y otra para N+1 y los intermedios interpolan?
			Ap[i]' = S[i] * Y[i]'		  							 
		END EXPAND_BLOCK
		EXPAND_BLOCK(i IN 1,N-1)   // Hay P de 1 a N, N+1 es igual a la 1 del siguiente componente
		(g[i+1]+0.5*dg[i+1])*Uc[i+1] + Apc[i+1]*P[i+1] = (g[i]+0.5*dg[i])*Uc[i]+Apc[i]*P[i]  -- Ec cant mov
		END EXPAND_BLOCK
		
---------------------------------------------------------------------- VOLUMEN DISCRETO (De 1 a N)--------------------------------------------------------------
			
		EXPAND_BLOCK (i IN 1,N)  -- 1 a N 
		   --Ab[i] = 0.5*(S[i]+S[i+1])*dx*Combustion[i]   -- Revisar
			Ab[i] = IF((Combustion[i]+Combustion[i+1])==0) 0 ELSE 0.5*(S[i]+S[i+1])*dx
			
			
			dg[i] = Rho_P*rp0[i]*Ab[i] 			-- Hay que usar r0 no rp0 - Usarla da mas elementos en la non linear box
			rp0[i] = (Factor_a*P[i]**Exponente)
			r0[i] = rp0[i]*eta[i]    
			
			Rho[i] = P[i] /(R_gas*T[i])
			SoundSpeed[i]= sqrt(gamma*R_gas*T[i])
			MACH[i]= Uc[i]/SoundSpeed[i]         
			
			Tt[i] = T[i]*(1+0.5*(gamma-1)*MACH[i]**2) 
			pt[i] = P[i]*(Tt[i]/T[i])**(gamma/(gamma-1))	
		
			Apc[i] = 0.5*(Ap[i]+Ap[i+1])
			Uc[i] = 0.5*(g[i]+g[i+1])/(Rho[i]*Apc[i])
		
		END EXPAND_BLOCK 
		EXPAND_BLOCK(i IN 1,N)  
			g[i+1]*Tt[i+1]=g[i]*Tt[i]+dg[i]*Tc   ---Ec energia
			
		END EXPAND_BLOCK
---------------------------------------- Calculos de combustible(De momento no da buenos resultados, repasarla) ----------------------------
	
	Ab_total = SUM(i IN 1,N; Ab[i] )
	V_total = SUM(i IN 1,N; Ab[i]*0.5*(Y[i]+Y[i+1]))
	Masa_total = V_total*Rho_P
	LoadFraction = V_total/(L*MATH.PI*(D/2)**2)


-------------------------------------------------- Combustión erosiva en volumenes (De 1 a N) ----------------------------------------------------------------
		
		
		IF (Comb_Erosiva == TRUE) INSERT
			EXPAND_BLOCK(i IN 1,N)  -- 1 a N
			g0_i[i] = (Rho[i]*Uc[i])/(Rho_P*rp0[i])  -- Cuidado con las velocidades aqui
			g_i[i] = g0_i[i]*(Re0[i]/1000)**-0.125
			Re0[i] = (Rho_P*rp0[i]*(S[i]/MATH.PI))/Viscosidad_camara  -- El Dh esta modificado segun paper para evitar sobrepresiones
			eta_temp[i] =	1 + 0.023*(g_i[i]**0.8-gth**0.8) 
			eta[i] = IF(eta_temp[i]<1) 1 ELSE 1 + 0.023*(g_i[i]**0.8-gth**0.8)
			END EXPAND_BLOCK 
		ELSE
			EXPAND_BLOCK(i IN 1,N)  
			g_i[i] = 0
			g0_i[i] = 0
			Re0[i] = 0
			eta_temp[i] = 1
			eta[i] = 1
			END EXPAND_BLOCK 
		END IF
		
	------------ Datos puerto ---------------------------------------------------------------------------------------------
		Coord [1]= entrada.Coord
		Coord[N+1]= salida.Coord
		
		g[1] = entrada.g
		g[N+1] = salida.g
		--(g[N+1])/(Rho[N]*Ap[N+1])= salida.U
		
		P[1]= entrada.Pout
		P[N]= salida.Pin
		Tt[1] = entrada.Tt
		Tt[N+1]= salida.Tt 
	
		Y[1] = entrada.Y
		S[1] = entrada.S
		Ap[1] = entrada.Ap
		
		Y[N+1] = salida.Y
		S[N+1] = salida.S
		Ap[N+1] = salida.Ap

		r0[1] = entrada.r0out
		r0[N] = salida.r0in
		
	   --Combustion[1] = entrada.Combustion
		Combustion[N+1] = salida.Combustion
		
		entrada.dgout = dg[1]
		salida.dgin= dg[N] 
		entrada.Ucout= Uc[1]
		salida.Ucin = Uc[N]
		entrada.Apcout = Apc[1]
		salida.Apcin = Apc[N]
		salida.gprev = g[N]
				
END COMPONENT


	
	

