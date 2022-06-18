/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: MChSolido3
 CREATION DATE: 09/06/2022
-----------------------------------------------------------------------------------------*/
USE LPRES
USE MATH

PORT SRB2 SINGLE
EQUAL							 REAL Coord				   UNITS u_m				RANGE ZERO,Inf			"Coordinate"	
EQUAL                    REAL g                 UNITS u_kg_s         RANGE ZERO,Inf       "Mass flow"
EQUAL							 REAL Pt				     UNITS u_Pa				RANGE ZERO,Inf		   "Total temperature"
EQUAL							 REAL Tt				     UNITS u_K				RANGE ZERO,Inf		   "Total temperature"


CONTINUOUS
END PORT
-------------------------------------------------------------------------------------------------------------------------------------------------------
COMPONENT Pared2
	PORTS
		OUT SRB2 salida
	DATA
		REAL Tc = 3330
	CONTINUOUS 
	salida.Tt = Tc		
	salida.g = 0
	salida.Coord = 0
END COMPONENT
-----------------------------------------------------------------------------------------------------------------------------------------------------
COMPONENT Conversor2
	PORTS
		IN SRB2 entrada
	DATA
		-- Datos gas y propelente --
		REAL c_star = 1500
		REAL R_gas = 307.93
		REAL gamma = 1.33
		--- Geometria tobera ---
		REAL D0 = 0.0254
		REAL A_th = 8.00E-03
	CONTINUOUS		
		1 = ((entrada.g*c_star)/(A_th*entrada.Pt))  -- Condicion de igualdad de gastos entre combustion y tobera 
		 ------ Erosion cte de la garganta -----
		 /*
		 rth' = alfa			-- Valor inicial rth = D0/2
		 At=MATH.PI*(rth)**2
		 */	
		------------- ECUACIONES TOBERA  ----------------
			
END COMPONENT
-------------------------------------------------------------------------------------------------------------------------------------------------------
COMPONENT MChSolido3(INTEGER N = 3, BOOLEAN Comb_Erosiva = TRUE )
	PORTS
		IN  SRB2 entrada
		OUT SRB2 salida
	DATA
			
			-- DATOS GENERALES COMBUSTIBLE Y PARAMETROS COMBUSTION
			REAL gamma = 1.33
			REAL Rho_P = 1750							 "Densidad del propelente"
			REAL M_m = 26.144    UNITS u_g 	 "Masa molar del propelente"
			REAL R_gas = 307.93						 "Constante R del gas"
			REAL Tc = 3300 							 "Temperatura de combustion adiabatica"
			REAL Cp = 1005
			-- Ley de vieille  
			REAL Exponente =	0.4					 "Exponente de Vieille"
			REAL Factor_a =	1.27955E-05			 "Factor de Vieille"	
			-- Combustion erosiva
			REAL Viscosidad_camara =	1.80E-05	 "Viscosidad de camara"
		   REAL gth =	35								 "Parametro de erosion"
			-- DATOS GEOMETRIA
			TABLE_1D S_y_1 = { {0,  0.2,0.200001}, -- Y
							  {0.075,0.2,0} } "Tabla S(y) de la interfaz de entrada" -- S 
		
			TABLE_1D S_y_N = { {0, 0.004,0.0040001,0.2}, -- Y
							  {0.185,0.185,0,0} }  "Tabla S(y) de la interfaz de salida"-- S	
							  
			REAL L = 0. UNITS u_m 		  "Longitud de la seccion"
			REAL D = 0.0635 UNITS u_m	  "Diametro carcasa cilindrica"
			REAL D0 = 0.0254 UNITS u_m   "Diametro inicial geometria cilindrica"
	DECLS
				------ Interfaces (1 a N+1)----------------

						-- Parametros geometricos
				DISCR REAL dx	
				REAL Y[N+1]
				REAL S[N+1]
				REAL Ap[N+1]
				REAL Coord [N+1]  UNITS "m"	
				REAL g[N+1]	
				REAL U[N+1]
				DISCR REAL Combustion[N+1]
				DISCR REAL A_carcasa
				---- "Volumenes" discretos	------ (1 a N)		
				REAL Ab[N] 				//Area de quemado del volumen
				REAL dg[N]				//Adicion de masa del volumen
				REAL P[N]        UNITS u_Pa				RANGE ZERO,Inf		   "Pressure"      
				REAL T[N]			UNITS u_Pa				RANGE ZERO,Inf		   "Temperature"  
				REAL Uc[N]
				REAL Apc[N]
		      -----------------------
				REAL rp0[N]				//Velocidad de recesion sin erosion
				REAL r0[N]				//Velocidad de recesion con erosion
				REAL Rho[N] 		
				REAL SoundSpeed[N]
				REAL MACH[N]
				REAL Tt[N]
				REAL Pt[N]
			   ------Calculos de combustible---------
					REAL Ab_total
					REAL V_total
					REAL Masa_total
					REAL LoadFraction
				-- Parametros de combustion erosiva --
				-- Modelo de mukunda --
					REAL g_i[N]
					REAL g0_i[N]  
					REAL Re0[N]
					REAL eta_temp[N]
					REAL eta[N]
								
	INIT		-----------Configuracion geometrica inicial, esto es provisional, se tendria que introducir en el experimento?-----------
		dx = L/N	
		FOR(i IN 1,N+1)
			Y[i] = 0
			--S[i] = 2*MATH.PI*0.5*D0				// De momento los valores iniciales se suponen los de un cilindro
		   Ap[i] = MATH.PI*(0.5*D0)**2
			A_carcasa = MATH.PI*(D/2)**2
			Combustion[i] = 1
		-- Aqui habria que implementar una distribucion de presiones inicial a partir de P0 = (rho_p*a*c_estrella* Abtotal/At)**(1/(1-n)) Quizas con una 
		-- variable intermedia que inicialice todo P[:] a P0
		END FOR
	DISCRETE
		EXPAND(i IN 1,N+1)
			WHEN (Y[i] >= (D-D0)/2) TOL 1E-05 THEN
				Combustion[i] = 0			// La combustion en una interfaz acaba cuando se llega a la carcasa
			END WHEN
			WHEN(SUM(i IN 1,N+1; Combustion[i]) <= N*0.6 ) THEN
		      STOP "******************NO FUEL LEFT*********************" // Realmente es que queden 3 interfaces o menos por quemar
				// Se podria hacer cuando la masa de propelente sea un 1-5% de la inicial (reg stacionario)?
			END WHEN
	CONTINUOUS		
		------------ Datos puerto ---------------------------------------------------------------------------------------------
		Coord [1]= entrada.Coord
		Coord[N+1]= salida.Coord
		g[1] = entrada.g
		g[N+1] = salida.g
		--U[1] = entrada.U
		--U[N+1] = salida.U
		
		Pt[1]= entrada.Pt
		Pt[N]= salida.Pt
		Tt[1] = entrada.Tt
		Tt[N]= salida.Tt 
		----------------------------------------------------------------INTERFACES (De 1 a N+1)---------------------------------------------------------------------
		
		EXPAND_BLOCK(i IN 1,N)
			Coord[i+1] = Coord[i] + dx			 
			g[i+1] = g[i]+dg[i]
			U[i] = g[i]/(Rho[i]*Ap[i]) 
		END EXPAND_BLOCK
			U[N+1] = g[N+1]/(Rho[N]*Ap[N+1])
		EXPAND_BLOCK(i IN 1,N-1)   
		(g[i+1]+0.5*dg[i+1])*Uc[i+1] + Apc[i+1]*P[i+1] = (g[i]+0.5*dg[i])*Uc[i]+Apc[i]*P[i]-- Ec cant mov
		END EXPAND_BLOCK
		------------------ Geometria -------------------
			Y[1]'	 = (2*r0[1]-Y[2]')*Combustion[1]				--- Modificar mediante una aproximacion mejor
			Y[N+1]'= (2*r0[N]-Y[N]')*Combustion[N+1]
		EXPAND_BLOCK(i IN 2,N)
			Y[i]' = 	0.5*(r0[i]+r0[i-1])*Combustion[i]
		END EXPAND_BLOCK
			S[1] = interpd1D(S_y_1,Y[1]) 
			S[N+1] = interpd1D(S_y_N,Y[N+1]) 
		EXPAND_BLOCK(i IN 2,N)
			S[i] = dx*(i-1)*(interpd1D(S_y_N,Y[i])-interpd1D(S_y_1,Y[i]) )/L + S[1]    
		END EXPAND_BLOCK
		EXPAND_BLOCK(i IN 1,N+1) 
			Ap[i]' = S[i] * Y[i]'
		END EXPAND_BLOCK

---------------------------------------------------------------------- VOLUMEN DISCRETO (De 1 a N)--------------------------------------------------------------
			
		EXPAND_BLOCK (i IN 1,N)  -- 1 a N 
			Ab[i] = IF((Combustion[i]+Combustion[i+1])==0) 0 ELSE 0.5*(S[i]+S[i+1])*dx
			
			dg[i] = Rho_P*rp0[i]*Ab[i] 			-- Cambiar rp0 por r0 
			rp0[i] = (Factor_a*P[i]**Exponente)
			r0[i] = rp0[i]*eta[i]    
			
			Rho[i] = P[i] /(R_gas*T[i])
			SoundSpeed[i]= sqrt(gamma*R_gas*T[i])
			MACH[i]= Uc[i]/SoundSpeed[i]
			
			Tt[i] = T[i]*(1+0.5*(gamma-1)*MACH[i]**2) 
			Pt[i] = P[i]*(Tt[i]/T[i])**(gamma/(gamma-1))	
		
			Apc[i] = 0.5*(Ap[i]+Ap[i+1])
			Uc[i] = 0.5*(g[i]+g[i+1])/(Rho[i]*Apc[i])
		END EXPAND_BLOCK 
		EXPAND_BLOCK(i IN 1,N-1)  
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
			eta[i] = IF(eta_temp[i]<1) 1 ELSE eta_temp[i]
			END EXPAND_BLOCK 
		ELSE
			EXPAND_BLOCK(i IN 1,N)  
			g_i[i] = 0
			g0_i[i] = 0
			Re0[i] = 0
			eta_temp[i] = 0
			eta[i] = 1
			END EXPAND_BLOCK 
		END IF		
END COMPONENT