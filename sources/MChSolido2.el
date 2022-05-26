/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: MChSolido2
 AUTHOR: Álvaro González Villarreal
 DESCRIptION: Seccion de motor cohete de propelente solido V0.2
 CREATION DATE: 02/04/2022
-----------------------------------------------------------------------------------------*/

USE MATH VERSION "3.2"

COMPONENT MChSolido2(INTEGER N = 3, BOOLEAN Comb_Erosiva = TRUE )
	PORTS
		IN  test entrada
		OUT test salida
	DATA
			REAL gamma =	1.33
			-- DATOS GENERALES COMBUSTIBLE Y PARAMETROS COMBUSTION
			REAL Rho_P = 1750	
			REAL R_gas = 307.93
			REAL Tc = 3300 
			-- Ley de vieille  
			REAL Exponente =	0.4	
			REAL Factor_a =	1.27955E-05	
			-- Combustion erosiva
			REAL Viscosidad_camara =	1.80E-05	
		   REAL gth =	35	
			-- DATOS GEOMETRIA
			REAL L = 5 UNITS "m" 			
			REAL D = 1 UNITS "m"
	DECLS
				------ Interfaces (1 a N+1)----------------
				-- Parametros geometricos
				REAL dx	
				REAL Y[N+1]
				REAL S[N+1]
				REAL Ap[N+1]
				REAL U[N+1] 
				REAL Coord [N+1]  UNITS "m"	
				REAL g[N+1]
				
				---- "Volumenes" discretos	------ (1 a N)		
				REAL Ab[N] 
				REAL dg[N]
				REAL P[N]             
				REAL T[N]
		      -----------------------
				REAL rp0[N]
				REAL Rho[N] 
				REAL SoundSpeed[N]
				REAL MACH[N]

				REAL Tt[N]
				REAL pt[N]
				-- Parametros de combustion erosiva --
				-- Modelo de mukunda --
				REAL g_i[N]
				REAL g0_i[N]  
				REAL Re0[N]
				REAL eta_temp[N]
				REAL eta[N]
			
	
	INIT							-----------Configuracion geometrica inicial-----------
		FOR(i IN 1,N+1)
			Y[i] = 0
			S[i] = 2*MATH.PI*0.5*0.0254
		   Ap[i] = MATH.PI*(0.5*0.0254)**2
		END FOR
	CONTINUOUS	
		-- Condiciones de contorno
		------------ Datos puerto -------------
		
		
		
		g[1] = entrada.g
		U[1]= entrada.U
		Coord [1]= entrada.Coord
		
		g[N+1] = salida.g
		U[N+1]= salida.U
		Coord[N+1]= salida.Coord
		
		 Y[1] = entrada.Y
		 S[1] = entrada.S
		Ap[1] = entrada.Ap
		--Y[1]' = entrada.Y'
		--S[1]' = entrada.S'
		--Ap[1]' = entrada.Ap'

		Y[N+1] = salida.Y
		S[N+1] = salida.S
		Ap[N+1] = salida.Ap
		--Y[N+1]' = salida.Y'
		--S[N+1]' = salida.S'
		--Ap[N+1]' = salida.Ap'
		
		
		
		P[1]= entrada.Pout
		T[1]= entrada.Tout
		eta[1]=entrada.etaout
		P[N]= salida.Pin
		T[N]= salida.Tin 
		eta[N]=salida.etain
		
		----------------------------------------------------------------INTERFACES (De 1 a N+1)---------------------------------------------------------------------
		dx = L/N
		EXPAND_BLOCK(i IN 1,N)
			Coord[i+1] = Coord[i] + dx			 
			g[i+1] = g[i]+dg[i]
			U[i+1] = (g[i+1])/(Rho[i]*Ap[i+1])   -- ¿Rho deberia ser la media entre volumenes adyacentes?
		END EXPAND_BLOCK	
		EXPAND_BLOCK(i IN 2,N)
			Y[i]' = 	0.5*(rp0[i]*eta[i]+rp0[i-1]*eta[i-1]) 
			S[i]'  =	2*MATH.PI*Y[i]'-- 0.5*(rp0[i]+rp0[i-1]) --Y[i]'  --*eta[i]
			Ap[i]' = S[i] * Y[i]'--    0.5*(rp0[i]+rp0[i-1])   --*eta[i]
		END EXPAND_BLOCK
---------------------------------------------------------------------- VOLUMEN DISCRETO (De 1 a N)--------------------------------------------------------------
		
		EXPAND_BLOCK (i IN 1,N)  -- 1 a N 
			Ab[i] = 0.5*(S[i]+S[i+1])*dx			
			dg[i] = Rho_P*rp0[i]*Ab[i] 			-- ¿Añadir terminos de masa de ignicion?
			rp0[i] = (Factor_a*P[i]**Exponente)
			
			Rho[i] = P[i] /(R_gas*T[i])
			SoundSpeed[i]= sqrt(gamma*R_gas*T[i])
			MACH[i]= (0.5*(U[i]+U[i+1]))/SoundSpeed[i]    -- Hacer velocidad media? Y la velocidad del sonido es la media?
			Tt[i] = T[i]*(1+0.5*(gamma-1)*MACH[i]**2) 
			pt[i] = P[i]*(Tt[i]/T[i])**(gamma/(gamma-1))
			
		END EXPAND_BLOCK 
		EXPAND_BLOCK(i IN 1,N-1)  -- 1 a N-1 define todos los volumenes mediante diferencias entre ellos, es necesaria una presion en un extremo?
			P[i+1] = P[i]  - ((g[i])*(U[i+1]-U[i]) + (dg[i])*U[i+1])/(0.5*(Ap[i]+Ap[i+1]))		--Algo falla en las presiones, un orden de mag por debajo de lo esperado
			--P[i+1] = P[i]-((Rho[i]*U[i]*Ap[i])*(U[i+1]-U[i])+g[i+1]*U[i+1])/(0.5*(Ap[i+1]+Ap[i]))
			
			Tt[i+1] = (g[i]*Tt[i]+dg[i]*Tc)/(g[i+1])	
			--T[i+1] = T[i] - ((0.5*(gamma-1))/(gamma*R_gas))  *  (U[i+1]**2 - U[i]**2) 
		END EXPAND_BLOCK
	
	
------------------------------------------------------------------------ Combustión erosiva en interfaces (De 1 a N+1) ----------------------------------------------------------------
		
		
		IF (Comb_Erosiva == TRUE) INSERT
			EXPAND_BLOCK(i IN 1,N)  -- 1 a N
			g0_i[i] = (Rho[i]*0.5*(U[i]+U[i+1]))/(Rho_P*rp0[i])
			g_i[i] = g0_i[i]*(Re0[i]/1000)**-0.125
			Re0[i] = (Rho_P*rp0[i]*(S[i]/MATH.PI))/Viscosidad_camara  -- El diametro hidraulico se puede definir con el volumen asi no se necesita S solo Ab, 4*Volumen/Ab
			eta_temp[i] =	1 + 0.023*(g_i[i]**0.8-gth**0.8) 
			eta[i] = IF(eta_temp[i]<1) 1 ELSE eta_temp[i]
			END EXPAND_BLOCK 
		ELSE
			EXPAND_BLOCK(i IN 1,N)  -- 1 a N
			g_i[i] = 0
			g0_i[i] = 0
			Re0[i] = 0
			eta_temp[i] = 0
			eta[i] = 1
			END EXPAND_BLOCK 
		END IF
		
	
	
	
	
	
	
	
	
	
			
END COMPONENT


	
	

