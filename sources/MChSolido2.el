/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: MChSolido2
 AUTHOR: Álvaro González Villarreal
 DESCRIptION: Seccion de motor cohete de propelente solido V0.2
 CREATION DATE: 02/04/2022
-----------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------
TAREAS A IMPLEMENTAR:
- ¿Establecer una distribucion inicial de P a partir del modelo cero dimensional? 

NOTAS

- Al hacer la partición es importante que la distribucion a introducir sea la de presiones

-----------------------------------------------------------------------------------------*/
USE MATH VERSION "3.2"

PORT test SINGLE
		
SUM	                   REAL g                 UNITS u_kg_s         RANGE ZERO,Inf       "Mass flow"
EQUAL		                REAL Ap                UNITS u_m2           RANGE ZERO,Inf       "Port Area"
EQUAL							 REAL S					   UNITS u_m				RANGE ZERO,Inf			"Perimeter of port"
EQUAL		                REAL U                 UNITS u_m2           RANGE ZERO,Inf       "Flow speed"			
EQUAL							 REAL Coord				   UNITS u_m				RANGE ZERO,Inf			"Coordinate"

EQUAL							 REAL Pout				   UNITS u_m				RANGE ZERO,Inf		   "Pre-interface pressure"
EQUAL							 REAL Pin				   UNITS u_m				RANGE ZERO,Inf		   "Post-interface pressure"						
-- Prueba


								 REAL rp0_intermedia						
CONTINUOUS
rp0_intermedia = (1.27955E-05)* (0.5*(Pin+Pout))**0.4    -- ¿Hay que pasar variables de la ley de vieille?
S'  =	2*MATH.PI*rp0_intermedia
Ap' = S*rp0_intermedia
END PORT

COMPONENT MChSolido2(INTEGER nodos = 3, BOOLEAN Comb_Erosiva = TRUE )
	PORTS
		IN  test  entrada
		OUT test salida
	DATA
			-- DATOS GENERALES COMBUSTIBLE Y PARAMETROS COMBUSTION
			REAL Rho_P = 1750	
			REAL gamma =	1.33
			REAL R_gas = 307.93
			REAL Tc = 3300 
			REAL c_star =	1498.7	
			-- Combustion erosiva
			REAL Viscosidad_camara =	1.80E-05	
			REAL gth =	35	
			-- Ley de vieille
			REAL Exponente =	0.4	
			REAL Factor_a =	1.27955E-05	
			-- DATOS GEOMETRIA
			REAL L = 5 UNITS "m" 			
	DECLS
				------ Interfaces (1 a N+1)----------------
				-- Parametros geometricos
				REAL dx	
				REAL S[nodos+1]
				REAL Ap[nodos+1]
				REAL U[nodos+1] 
				REAL Coord [nodos+1]  UNITS "m"	
				REAL g[nodos+1]
				-- Parametros de combustion erosiva
				/* REAL Kp[nodos+1]
				REAL Jx[nodos+1]
				REAL g0[nodos+1]  
				REAL Re[nodos+1]
				REAL eta_temp[nodos+1]
				REAL eta[nodos+1] */
			------ "Volumenes" discretos	------ (1 a N)		
				-- Parametros termodinamicos
				REAL Ab[nodos] 
				REAL dg[nodos]
				REAL P[nodos]             
				REAL T[nodos]
		      -----------------------
				REAL rp0[nodos]
				REAL Rho[nodos] 
				REAL SoundSpeed[nodos]
				REAL MACH[nodos]
				-- Parametros de remanso
				REAL Tt[nodos]
				REAL pt[nodos]
				
	
	INIT							-----------Configuracion geometrica inicial-----------
		FOR(i IN 1,nodos+1)
			S[i] = 4.09E-01
			Ap[i] = 1.33E-02
		END FOR
	CONTINUOUS	
		-- Condiciones de contorno
		------------ Datos puerto (WIP)-------------
		g[1] = entrada.g
		Ap[1] = entrada.Ap
		S[1] = entrada.S
		U[1]= entrada.U
		Coord [1]= entrada.Coord
		
		g[nodos+1] = salida.g
		Ap[nodos+1] = salida.Ap
		S[nodos+1] = salida.S
		U[nodos+1]= salida.U
		Coord[nodos+1]= salida.Coord
		
		P[1]= entrada.Pout
		P[nodos]= salida.Pin 
		
		---------------------------------------------
		
		
		-- Condicion de igualdad de gastos entre combustion y tobera --  ¿Tramos intermedios? 
		--1 = (g[nodos+1]*c_star)/(A_th*pt[nodos])  Habria que imponerla en los boundary (presion de remanso en el inicio)-- 
		
		
---------------------------------------------------------------------------INTERFACES (De 1 a N+1)---------------------------------------------------------------------
		dx = L/(nodos)
		EXPAND_BLOCK(i IN 2,nodos+1)
			Coord[i] = Coord[i-1] + dx			 
			g[i] = g[i-1]+dg[i-1]
			U[i] = (g[i])/(Rho[i-1]*Ap[i])   -- ¿Rho deberia ser la media entre volumenes adyacentes?
		END EXPAND_BLOCK	
		EXPAND_BLOCK(i IN 2,nodos)
			S[i]'  =	2*MATH.PI*   0.5*(rp0[i]+rp0[i-1])   --*eta[i]
			Ap[i]' = S[i]*        0.5*(rp0[i]+rp0[i-1])   --*eta[i]
		END EXPAND_BLOCK
---------------------------------------------------------------------- VOLUMEN DISCRETO (De 1 a N)--------------------------------------------------------------
		
		EXPAND_BLOCK (i IN 1,nodos)  -- 1 a N 
		Ab[i] = 0.5*(S[i]+S[i+1])*dx
		dg[i] = Rho_P*rp0[i]*Ab[i] 			-- ¿Añadir terminos de masa de ignicion?
		rp0[i] = Factor_a*P[i]**Exponente            
		END EXPAND_BLOCK 
		EXPAND_BLOCK(i IN 1,nodos-1)  
		P[i+1] = P[i]-((g[i])*(U[i+1]-U[i]) + dg[i]*U[i+1])/Ap[i+1]		
		Tt[i+1] = (g[i]*Tt[i]+dg[i]*Tc)/(g[i+1])	
		END EXPAND_BLOCK
		EXPAND_BLOCK(i IN 1,nodos) -- Variables de remanso -- 
		T[i] = Tt[i]*(1+0.5*(gamma-1)*MACH[i]**2)**(-1)
		P[i] = pt[i]*(1+0.5*(gamma-1)*MACH[i]**2)**(-(gamma/(gamma-1))) --P[i]*(Tt[i]/T[i])**(gamma/(gamma-1)) 
		Rho[i] = P[i] /(R_gas*T[i])
		SoundSpeed[i]= sqrt(gamma*R_gas*T[i])
		MACH[i]= 0.5*(U[i]+U[i+1])/SoundSpeed[i]    -- Hacer velocidad media? Y la velocidad del sonido es la media?
		
		END EXPAND_BLOCK			
	
------------------------------------------------------------------------ Combustión erosiva (De 1 a N+1) ----------------------------------------------------------------
		
	/*	
		IF (Comb_Erosiva == TRUE) INSERT
			EXPAND_BLOCK(i IN 1,nodos)  -- 1 a N
			Kp[i] = Ap[i]/Ab[i]
			Jx[i] = Kp[i]/Klemmug
			g0[i] = (g[i]/Ap[i])/(Rho_P*rp0[i])
			Re[i] = ((Rho_P*rp0[i]*  4*Ap[i])/S[i]  )/Viscosidad_camara  -- El diametro hidraulico se puede definir con el volumen asi no se necesita S solo Ab, 4*Volumen/Ab
			eta_temp[i] =	1 + 0.023*(g0[i]**0.8-gth**0.8) 
			eta[i] = IF(eta_temp[i]<1) 1 ELSE eta_temp[i]
			END EXPAND_BLOCK 
			Kp[nodos+1] = Ap[nodos+1]/Ab[nodos]  -- La ultima interfaz toma los datos del ultimo volumen ¿
			Jx[nodos+1] = Kp[nodos+1]/Klemmug
			g0[nodos+1] = (g[nodos+1]/Ap[nodos+1])/(Rho_P*rp0[nodos+1])
			Re[nodos+1] = ((Rho_P*rp0[nodos+1]*4*Ap[nodos+1])/S[nodos+1])/Viscosidad_camara
			eta_temp[nodos+1] =	1 + 0.023*(g0[nodos+1]**0.8-gth**0.8) 
			eta[nodos+1] = IF(eta_temp[nodos+1]<1) 1 ELSE eta_temp[nodos+1]
		ELSE
			EXPAND_BLOCK(i IN 1,nodos+1)  -- 1 a N
			Kp[i] = 0
			Jx[i] = 0
			g0[i] = 0
			Re[i] = 0
			eta_temp[i] = 0
			eta[i] = 1
			END EXPAND_BLOCK 
		END IF
		
	*/
	
	
	
	
	
	
	
	
			
END COMPONENT


	
	

