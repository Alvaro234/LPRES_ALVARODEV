/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: Pruebageometria
 CREATION DATE: 29/05/2022
-----------------------------------------------------------------------------------------*/
USE MATH

COMPONENT SimuladorGeometria(INTEGER N = 10)
	DATA
		
		TABLE_1D S1 = {{0,  0.02}, -- Y
							  {0.075,0.2} }  -- S
		
		TABLE_1D S2 = { {0, 0.004,0.004000000000000001,0.004000000000000002}, -- Y
							  {0.185,0.186,0,0} }  -- S					 
	
			-- Se interpolan los extremos del componente para tener 
		REAL L = 0.708 UNITS u_m 		  "Longitud de la seccion"
		REAL D = 0.0635 UNITS u_m	  "Diametro carcasa cilindrica"
		REAL D0 = 0.0240 UNITS u_m   "Diametro inicial geometria cilindrica"
		REAL Rho_P = 1750
	
	
	DECLS
		DISCR REAL dx
		REAL Y[N+1]
	   REAL S[N+1]
		REAL Ap[N+1]
		REAL Ab[N]
		DISCR REAL Combustion[N+1]
		REAL Ab_total
		REAL V_total
		REAL Masa_total
		REAL LoadFraction
	INIT		
			dx = L/N	
			FOR(i IN 1,N+1)
				Y[i] = 0
				Ap[i] = MATH.PI*(0.5*D0)**2   -- No es del todo correcto, cada nodo tendra un area en funcion de su geometria
				Combustion[i] = 1	
			END FOR
	DISCRETE
			EXPAND(i IN 1,N+1)
				WHEN (Y[i] >= (D-D0)/2) TOL 1E-05 THEN
					Combustion[i] = 0			// La combustion en una interfaz acaba cuando se llega a la carcasa
				END WHEN
				WHEN(SUM(i IN 1,N+1; Combustion[i]) <= N*0.1 ) THEN
			      STOP "******************NO FUEL LEFT*********************" // Realmente es que queden 3 interfaces o menos por quemar
					// Se podria hacer cuando la masa de propelente sea un 1-5% de la inicial (reg stacionario)?
				END WHEN	
				
				
	CONTINUOUS 
		S[1] = interpd1D(S1,Y[1]) 
		S[N+1] = interpd1D(S2,Y[N+1])
		EXPAND_BLOCK(i IN 2,N)
		S[i] = dx*(i-1)*(interpd1D(S2,Y[i])-interpd1D(S1,Y[i]))/L + S[1]    
		END EXPAND_BLOCK
		 
		
		EXPAND_BLOCK(i IN 1,N+1)
		Y[i]' = Combustion[i]*(0.00001*i**2+0.01)
		Ap[i]' = S[i] * Y[i]'
		END EXPAND_BLOCK	
		
		
		EXPAND_BLOCK(i IN 1,N)
		Ab[i] = IF((Combustion[i]+Combustion[i+1])==0) 0 ELSE 0.5*(S[i]+S[i+1])*dx   
		END EXPAND_BLOCK
		
		Ab_total = SUM(i IN 1,N; Ab[i] )
		V_total = SUM(i IN 1,N; Ab[i]*0.5*(Y[i]+Y[i+1]))
		Masa_total = V_total*Rho_P
		LoadFraction = V_total/(L*MATH.PI*(D/2)**2)

END COMPONENT
