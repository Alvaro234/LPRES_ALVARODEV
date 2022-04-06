/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: MChSolido2
 AUTHOR: Álvaro González Villarreal
 DESCRIptION: Seccion de motor cohete de propelente solido V0.2
 CREATION DATE: 02/04/2022
-----------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------
TAREAS A IMPLEMENTAR:
- Establecer una distribucion inicial de P a partir del modelo cero dimensional? 
- Caida de presion con coeficiente landa entre el nodo N-1 (final del grano) hasta la tobera N ¿?

NOTAS

- Al hacer la partición es importante que la distribucion a introducir sea la de presiones y sin mezclar con la velocidad

-----------------------------------------------------------------------------------------*/
USE MATH VERSION "3.2"

PORT test SINGLE
		
		             REAL T                 UNITS u_K            RANGE ZERO,Inf       "Total temperature"
	                REAL P                 UNITS u_Pa           RANGE ZERO,Inf       "Total pressure"
		             REAL g                 UNITS u_kg_s         RANGE ZERO,Inf       "Mass flow"
	                REAL Ap              UNITS u_m2           RANGE ZERO,Inf       "Port Area"
						 REAL S					UNITS u_m			"Perimeter of port"
--EQUAL	             REAL c_star               UNITS u_m_s          RANGE ZERO,Inf       "Characteristic velocity"
--EQUAL			       REAL A_th                 UNITS u_m2           RANGE ZERO,Inf       "Throat area"

END PORT

COMPONENT MChSolido2(INTEGER nodos = 3, BOOLEAN Comb_Erosiva = TRUE )
	PORTS
		IN test  entrada
		OUT test salida
	DATA
			-- DATOS GENERALES COMBUSTIBLE Y PARAMETROS COMBUSTION
			REAL Rho_P = 1750	
			REAL gamma =	1.33	
			REAL Viscosidad_camara =	1.80E-05	
			REAL gth =	35	
			REAL Temperatura_de_combustion =3300 
			REAL c_star =	1498.7	
			REAL Exponente =	0.4	
			REAL Factor_a =	1.27955E-05	
			REAL R_gas = 307.93  		
			-- DATOS MOTOR
			REAL L = 5.48E+00 UNITS "m" 
			-- DATOS TOBERA 
			REAL A_th = 8.00E-03	UNITS "m2"
			REAL Klemmug = 280.64				
	DECLS
			CONST REAL Pi = MATH.PI    -- eliminar, ya esta en MATH
			
			------ Interfaces (1 a N+1)----------------
				-- Parametros geometricos
				REAL dx	
				REAL S[nodos+1]
				REAL Ap[nodos+1]
				REAL rp0[nodos+1]
				REAL Coord [nodos+1]  UNITS "m"	
				REAL g[nodos+1]
				-- Parametros de combustion erosiva
				REAL Kp[nodos+1]
				REAL Jx[nodos+1]
				REAL g0[nodos+1]  
				REAL Re[nodos+1]
				REAL eta_temp[nodos+1]
				REAL eta[nodos+1]
			------ "Volumenes" discretos	------ (1 a N)		
				-- Parametros termodinamicos
				BOUND REAL P[nodos]             
				REAL U[nodos] 
				REAL T[nodos]
		      -----------------------
				REAL Rho[nodos] 
				REAL SoundSpeed[nodos]
				REAL MACH[nodos]
				-- Parametros de remanso
				REAL Tt[nodos]
				REAL pt[nodos]
				REAL Ab[nodos] 
				REAL dg[nodos]
	
	INIT							-----------Configuracion geometrica inicial-----------
		FOR(i IN 1,nodos+1)
			S[i] = 4.09E-01
			Ap[i] = 1.33E-02
		END FOR
		
	CONTINUOUS	
		-- Condiciones de contorno
		dx = L/(nodos-1)
		------------ Datos puerto (WIP)-------------
		g[1] = 0
		U[1]= 0
		T[1] = Temperatura_de_combustion
		--P[1] = (Rho_P*Factor_a*c_star*(Ab[1]/A_th))    
		
		T[nodos] = salida.T
		P[nodos] = salida.P
		g[nodos] = salida.g
		Ap[nodos+1] = salida.Ap
		S[nodos+1] = salida.S
		
		-- c_star = salida.c_star             
      -- A_th = salida.A_th
		
		T[1] = entrada.T
		--P[1] = entrada.P
		g[1] = entrada.g
		Ap[1] = entrada.Ap
		S[1] = entrada.S
		
		
		
		-- Condicion de igualdad de gastos entre combustion y tobera --
		1 = (g[nodos]*c_star)/(A_th*pt[nodos])    
		
		
		----------------------------------------------------------------------INTERFACES (De 1 a N+1)---------------------------------------------------------------------
		
		-- Geometria de las interfaces en el tiempo -- Combustión erosiva
		-- Posicion de las interfaces y adicion de gasto por combustion en volumen previo
		Coord[1]=0
		EXPAND_BLOCK(i IN 2,nodos+1)
		Coord[i] = Coord[i-1] + dx	
		g[i] = g[i-1]+dg[i-1]
		END EXPAND_BLOCK	
		IF(Comb_Erosiva == TRUE) INSERT
			EXPAND_BLOCK(i IN 1,nodos)
				S[i]'  =	2*Pi*rp0[i]*eta[i]
				Ap[i]' = S[i]*rp0[i]*eta[i]
				rp0[i] = Factor_a*P[i]**Exponente
			END EXPAND_BLOCK
				S[nodos+1]'  =	2*Pi*rp0[nodos+1]*eta[nodos+1]
				Ap[nodos+1]' = S[nodos+1]*rp0[nodos+1]*eta[nodos+1]
				rp0[nodos+1] = Factor_a*P[nodos]**Exponente
		ELSE
			EXPAND_BLOCK(i IN 1,nodos)
				S[i]'  =	2*Pi*rp0[i]
				Ap[i]' =  S[i]*rp0[i]
				rp0[i] = Factor_a*(P[i])**Exponente
			END EXPAND_BLOCK
				S[nodos+1]'  =	2*Pi*rp0[nodos+1]
				Ap[nodos+1]' = S[nodos+1]*rp0[nodos+1]
				rp0[nodos+1] = Factor_a*P[nodos]**Exponente
		END IF
		
		
		----------------------------------------------------------------- VOLUMEN DISCRETO (De 1 a N)--------------------------------------------------------------
		
		EXPAND_BLOCK (i IN 1,nodos)  -- 2 a N 
		Ab[i] = 0.5*(S[i]+S[i+1])*dx
		dg[i] = 0.5*Rho_P*( rp0[i]*S[i] + rp0[i+1]*S[i+1] )*dx  --Añadir terminos de masa de ignicion?
		END EXPAND_BLOCK 
		EXPAND_BLOCK(i IN 1,nodos-1)   -- El nodo de entrada dependera de la conexion en el puerto, si no hay ninguna --->> P[1] = (Rho_P*Factor_a*c_star*(Ab[1]/A_th))
		P[i+1] = P[i]-((g[i])*(U[i+1]-U[i])+g[i+1]*U[i+1])/(0.5*(Ap[i+1]+Ap[i]))
		U[i+1] = (g[i+1])/(Rho[i+1]*Ap[i+1])	
		T[i+1] = T[i] - ((0.5*(gamma-1))/(gamma*R_gas))  *  (U[i+1]**2 - U[i]**2) 
		END EXPAND_BLOCK
		EXPAND_BLOCK(i IN 1,nodos) -- Variables de remanso
		Tt[i] = T[i]*(1+0.5*(gamma-1)*MACH[i]**2)
		pt[i] = P[i]*(Tt[i]/T[i])**(gamma/(gamma-1)) 
		SoundSpeed[i]= sqrt(gamma*R_gas*T[i])
		MACH[i]= U[i]/SoundSpeed[i]
		Rho[i] = P[i] /(R_gas*T[i])
		END EXPAND_BLOCK			
	
		------------------------------------------------------------------ Combustión erosiva  ----------------------------------------------------------------
		IF (Comb_Erosiva == TRUE) INSERT
			EXPAND_BLOCK(i IN 1,nodos)  -- 1 a N
			Kp[i] = Ap[i]/Ab[i]
			Jx[i] = Kp[i]/Klemmug
			g0[i] = (g[i]/Ap[i])/(Rho_P*rp0[i])
			Re[i] = ((Rho_P*rp0[i]*4*Ap[i])/S[i])/Viscosidad_camara
			eta_temp[i] =	1 + 0.023*(g0[i]**0.8-gth**0.8) 
			eta[i] = IF(eta_temp[i]<1) 1 ELSE eta_temp[i]
			END EXPAND_BLOCK 
			Kp[nodos+1] = Ap[nodos+1]/Ab[nodos]  -- La ultima interfaz toma los datos del ultimo volumen
			Jx[nodos+1] = Kp[nodos+1]/Klemmug
			g0[nodos+1] = (g[nodos+1]/Ap[nodos+1])/(Rho_P*rp0[nodos+1])
			Re[nodos+1] = ((Rho_P*rp0[nodos+1]*4*Ap[nodos+1])/S[nodos+1])/Viscosidad_camara
			eta_temp[nodos+1] =	1 + 0.023*(g0[nodos+1]**0.8-gth**0.8) 
			eta[nodos+1] = IF(eta_temp[nodos+1]<1) 1 ELSE eta_temp[nodos+1]
		ELSE
		END IF
		
			
END COMPONENT




	
	

