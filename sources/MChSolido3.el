/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: MChSolido3
 CREATION DATE: 09/06/2022
-----------------------------------------------------------------------------------------*/
USE LPRES
USE MATH

PORT SRB2 SINGLE
EQUAL							 REAL Coord				  UNITS u_m				RANGE ZERO,Inf			"Coordinate"	
EQUAL                    REAL g                UNITS u_kg_s       RANGE ZERO,Inf       "Mass flow"
EQUAL							 REAL Pt				     UNITS u_Pa			RANGE ZERO,Inf		   "Total temperature"
EQUAL							 REAL Tt				     UNITS u_K				RANGE ZERO,Inf		   "Total temperature"

END PORT

-------------------------------------------------------------------------------------------------------------------------------------------------------

COMPONENT Pared2
	PORTS
		OUT SRB2 salida
	DATA
		REAL Tc = 3330	 UNITS u_K	RANGE ZERO,Inf "Adiabatic combustion temperature"
	CONTINUOUS 
		-------------- Head End Boundary Conditions -----------------
		salida.Tt = Tc				
		salida.g = 0
		salida.Coord = 0
		-------------------------------------------------------------
END COMPONENT

-----------------------------------------------------------------------------------------------------------------------------------------------------
COMPONENT Conversor2
	PORTS
		IN SRB2 entrada
	DATA
		-- Datos gas y propelente --
		REAL c_star = 1500			UNITS u_m_s		"Velocidad caracteristica"
		REAL R_gas = 307.93			UNITS u_J_kgK	 "Constante R del gas"
		REAL gamma = 1.33				UNITS no_units	 "Relacion de calores especificos"
		--- Geometria tobera ---
		REAL D0 = 0.0254		UNITS u_m	"Diametro inicial de la tobera"
		REAL ER = 10			UNITS no_units "Relacion de areas"
		REAL Pa = 101325		UNITS u_Pa			   "Ambient pressure"
		---- Erosion garganta ---
		REAL alfa = 0.0002	UNITS u_m_s	"Velocidad de recesion de la garganta"
	DECLS
		------------- 
		REAL rth			"Recesion de la tobera"
		REAL A_th = 8.00E-03
		------------- Actuacion tobera ---------------
		DISCR REAL Fgamma			UNITS no_units
		REAL Ct						UNITS no_units
		REAL T						UNITS u_N
		REAL Pexit					UNITS u_Pa
		REAL TotalImpulse			UNITS u_N
		REAL ISP						UNITS u_s
	INIT
		Fgamma = sqrt(gamma)* ( 2/(gamma+1) )**( (gamma+1)/(2*(gamma-1)) )
		TotalImpulse = 0
		rth = D0/2
		A_th = MATH.PI*(D0/2)**2
	CONTINUOUS		
		c_star = (A_th*entrada.Pt)/entrada.g -- Condicion de igualdad de gastos entre combustion y tobera --
		------ Erosion lineal de la garganta -----
		 rth' = alfa			-- Valor inicial rth = D0/2
		 A_th=MATH.PI*(rth)**2
		------------- ECUACIONES TOBERA (Adaptada) ----------------    
		Pexit = Pa			// Tobera daptada
		Ct = Fgamma*sqrt((   2*gamma/(gamma-1) )*  (1- (Pexit/entrada.Pt)**((gamma-1)/gamma) )       )
		T = Ct*entrada.Pt*A_th
		TotalImpulse' = T
		ISP = TotalImpulse/(entrada.g*9.81)
		--ER = Fgamma/ (( (Pexit/entrada.Pt)**(1/gamma)) * sqrt((   2*gamma/(gamma-1) )*  (1- (Pexit/entrada.Pt)**((gamma-1)/gamma) )  ))	
	
END COMPONENT
-------------------------------------------------------------------------------------------------------------------------------------------------------
COMPONENT MChSolido3(INTEGER N = 3, BOOLEAN Comb_Erosiva = TRUE )
	PORTS
		IN  SRB2 entrada
		OUT SRB2 salida
	DATA
			-- DATOS GENERALES COMBUSTIBLE Y PARAMETROS COMBUSTION
			REAL gamma = 1.33		UNITS no_units	 "Relacion de calores especificos"
			REAL Rho_P = 1750		UNITS u_kg_m3	 "Densidad del propelente"
			REAL M_m = 26.144    UNITS u_g 	 	 "Masa molar del propelente"
			REAL R_gas = 307.93	UNITS u_J_kgK	 "Constante R del gas"
			REAL Tc = 3300 		UNITS u_K		 "Temperatura de combustion adiabatica"
			-- Ley de vieille  
			REAL Exponente =	0.4					 "Vieille law pressure exponent"
			REAL Factor_a =	1.27955E-05			 "Vieille law coefficient"	
			-- Combustion erosiva
			REAL Viscosidad_camara =	1.80E-05	 "Viscosidad de camara"
		   REAL gth =	35								 "Parametro de erosion"
			-- DATOS GEOMETRIA
			TABLE_1D S_y_1 = { {0, 2}, 
							  		 {1, 2} } UNITS u_m	"Tabla S(y) de la interfaz de entrada"
		
			TABLE_1D S_y_N = { {0, 2}, 
							  		 {1, 2} } UNITS u_m 	"Tabla S(y) de la interfaz de salida"
							  
			REAL L = 0. 			UNITS u_m 		  "Longitud de la seccion"
			REAL D = 0.0635 		UNITS u_m	     "Diametro carcasa cilindrica"
			REAL D0 = 0.0254 		UNITS u_m        "Diametro inicial geometria cilindrica"
	DECLS
				------ Interfaces (1 a N+1)----------------

						-- Parametros geometricos
				DISCR REAL dx		UNITS u_m
				REAL Y[N+1]			UNITS u_m
				REAL S[N+1]			UNITS u_m
				REAL Ap[N+1]		UNITS u_m2
				REAL Coord [N+1]  UNITS u_m	
				REAL g[N+1]			UNITS u_kg_s
				REAL U[N+1]			UNITS u_m_s
				DISCR REAL Combustion[N+1]	UNITS no_units
				DISCR REAL A_carcasa			UNITS u_m2
				---- "Volumenes" discretos	------ (1 a N)		
				REAL Ab[N] 		  UNITS u_m2	//Area de quemado del volumen
				REAL dg[N]		  UNITS u_kg_s		//Adicion de masa del volumen
				REAL P[N]        UNITS u_Pa				RANGE ZERO,Inf		   "Pressure"      
				REAL T[N]		  UNITS u_K				RANGE ZERO,Inf		   "Temperature"  
				REAL Uc[N]		  UNITS u_m_s
				REAL Apc[N]		  UNITS u_m2
		      -----------------------				
				REAL rp0[N]		  UNITS u_m_s	//Velocidad de recesion sin erosion
				REAL r0[N]		  UNITS u_m_s	//Velocidad de recesion con erosion
				REAL Rho[N] 	  UNITS u_kg_m3
				REAL SoundSpeed[N] UNITS u_m_s
				REAL MACH[N]	  UNITS no_units
				REAL Tt[N]		  UNITS u_K
				REAL Pt[N]		  UNITS u_K
			    ------Calculos de combustible---------
					REAL Ab_total		UNITS u_m2
					REAL Vburnt			UNITS u_m3
					REAL Masaburnt		UNITS u_kg
					REAL LoadFraction UNITS no_units
				-- Parametros de combustion erosiva --
				-- Modelo de mukunda --
					REAL g_i[N]
					REAL g0_i[N]  
					REAL Re0[N]
					REAL eta_temp[N]
					REAL eta[N]
								
	INIT		
	   dx = L/N	   -- Malla equiespaciada --- 
		Vburnt = 0
		A_carcasa = MATH.PI*(D/2)**2
			FOR(i IN 1,N+1)
				Y[i] = 0
				Ap[i] = MATH.PI*(0.5*D0)**2 // De momento los valores de areas iniciales se suponen los de un cilindro
				Combustion[i] = 1
			END FOR
	DISCRETE
		EXPAND(i IN 1,N+1)
			WHEN (Y[i] >= (D-D0)/2) TOL 1E-06 THEN
				Combustion[i] = 0			// La combustion en una interfaz acaba cuando se llega a la carcasa
			END WHEN
			WHEN(SUM(i IN 1,N+1; Combustion[i]) <= N*0.3 ) THEN			-- Cuando la derivada de la masa sea negativa?
		      --STOP "******************NO FUEL LEFT*********************" // Realmente es que queden 3 interfaces o menos por quemar
			END WHEN
	CONTINUOUS		
		------------ Datos puerto ---------------------------------------------------------------------------------------------
		Coord [1]= entrada.Coord
		Coord[N+1]= salida.Coord
		g[1] = entrada.g
		g[N+1] = salida.g
		Pt[1]= entrada.Pt
		Pt[N]= salida.Pt
		Tt[1] = entrada.Tt
		Tt[N]= salida.Tt 
		----------------------------------------------------------------INTERFACES (De 1 a N+1)---------------------------------------------------------------------
		------------------ Geometria -------------------
			Y[1]'	 = (2*r0[1]-Y[2]')*Combustion[1]		-- Los extremos se aproximan mediante la velocidad de recesion del volumen y la interfaz anterior
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
		EXPAND_BLOCK(i IN 1,N)
			Coord[i+1] = Coord[i] + dx			 
			g[i+1] = g[i]+dg[i]
			U[i] = g[i]/(Rho[i]*Ap[i]) 
		END EXPAND_BLOCK
			U[N+1] = g[N+1]/(Rho[N]*Ap[N+1])
		EXPAND_BLOCK(i IN 1,N-1)   
		(g[i+1]+0.5*dg[i+1])*Uc[i+1] + Apc[i+1]*P[i+1] = (g[i]+0.5*dg[i])*Uc[i]+Apc[i]*P[i]-- Ec cant mov
		END EXPAND_BLOCK
		
---------------------------------------------------------------------- VOLUMEN DISCRETO (De 1 a N)--------------------------------------------------------------
		EXPAND_BLOCK (i IN 1,N)  -- 1 a N 
			Ab[i] = IF((Combustion[i]+Combustion[i+1])==0) 1e-6 ELSE 0.5*(S[i]+S[i+1])*dx
			
			dg[i] = Rho_P*r0[i]*Ab[i] 					 
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
		--T[i+1]*(g[i+1]/g[i]) = T[i]+Tc*(dg[i]/g[i]) -- Intento de correccion de fallo por final de combustion de la interfaz 1 y 2
		END EXPAND_BLOCK

------------------------------- Calculos de combustible(De momento no da buenos resultados (quizas por la falta de transitorio), repasarla) ----------
	
	Ab_total = SUM(i IN 1,N; Ab[i] )
	Vburnt' = SUM(i IN 1,N; Ab[i]*0.5*(Y[i]+Y[i+1]))
	Masaburnt = Vburnt*Rho_P
	LoadFraction = Vburnt/(L*MATH.PI*(D/2)**2)

-------------------------------------------------- Combustión erosiva Mukunda en volumenes (De 1 a N) ----------------------------------------------------------------
		IF (Comb_Erosiva == TRUE) INSERT
			EXPAND_BLOCK(i IN 1,N)  -- 1 a N
			g0_i[i] = (Rho[i]*Uc[i])/(Rho_P*rp0[i])  -- Cuidado con las velocidades aqui
			g_i[i] = g0_i[i]*(Re0[i]/1000)**-0.125
			Re0[i] = (Rho_P*rp0[i]*(0.5*(S[i]+S[i+1])/MATH.PI))/Viscosidad_camara  -- El Dh esta modificado segun paper para evitar sobrepresiones
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
------------------------------------------------- Combustiónes erosivas extras en volumenes (De 1 a N) ----------------------------------------------------------------
END COMPONENT
