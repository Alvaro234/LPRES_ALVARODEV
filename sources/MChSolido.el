/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: MChSolido
 AUTHOR: Álvaro González Villarreal
 DESCRIPTION: Seccion de motor cohete de propelente solido V0.1
 CREATION DATE: 19/03/2022
-----------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------
TAREAS A IMPLEMENTAR:

- Implementar evolución de Ab, Ap y S con el tiempo.
- Corregir funcionamiento de eta (Algoritmo de mukunda, coeficiente de erosion)/ ¿Implementarlo como función? eta >= 1 DONE
- Establecer una distribucion inicial a partir del modelo cero dimensional? 

- Quitar datos sobrantes del codigo

NOTAS

- Al hacer la partición es importante que la distribucion a introducir sea la de presiones y sin mezclar con la velocidad

-----------------------------------------------------------------------------------------*/
COMPONENT MChSolido(INTEGER nodos = 11)
	-- De momento sin puertos	
	DATA
			REAL Rho_P = 1750	UNITS "kg/m3"
			REAL gamma =	1.33	
			REAL Viscosidad_camara =	1.80E-05	
			BOOLEAN Combustion_erosiva = TRUE	
			REAL gth =	35	
			REAL Temperatura_de_combustion =3300 UNITS	"K"
			REAL ce_estrella =	1498.7	UNITS "m/s"
			REAL Exponente =	0.4	
			REAL Factor_a =	1.27955E-05	
			REAL R_gas = 307.93  UNITS "J/kgK"			
			-- DATOS MOTOR
			REAL Area_de_garganta	 = 8.00E-03	UNITS "m2"
			REAL Klemmug = 280.64	
			REAL Par_func_Y =	0.404
			REAL L = 5.48E+00 UNITS "m" 				
	DECLS
			CONST REAL Pi = 3.1415926
			REAL Coord [nodos]  UNITS "m"	
			REAL dx			
			REAL Ab[nodos] UNITS "m2"
			REAL Kp[nodos]
			REAL Jx[nodos]
			REAL g0[nodos]  UNITS "no_units"
			REAL gRe[nodos] UNITS "no_units"
			REAL Re[nodos]
			REAL eta[nodos]
			REAL g[nodos]
			REAL rp0[nodos] 
			REAL dg[nodos]
			REAL P[nodos] UNITS "Pa"
			REAL U[nodos] UNITS "m/s"
			REAL T[nodos] UNITS "K"
			REAL Rho[nodos] UNITS "kg/m3"
			REAL SoundSpeed[nodos]
			REAL MACH[nodos]
			REAL Tt[nodos]
			REAL Pt[nodos]
			REAL S[nodos] 
			REAL Ap[nodos] 
			REAL eta_temp[nodos]
	INIT
		FOR(i IN 1,nodos)
			S[i] = 4.09E-01
			Ap[i] = 1.33E-02
		END FOR
	CONTINUOUS	
		-- Condiciones de contorno
		dx = L/(nodos-1)
		Coord[1]=0
		Ab[1] = 0
		dg[1] = 0
		g[1] = 0
		U[1]= 0
		T[1] = Temperatura_de_combustion
		1 = (g[nodos]*ce_estrella)/(Area_de_garganta*Pt[nodos]) 
	
		
		EXPAND_BLOCK (i IN 1,nodos)
			Kp[i] = Ap[i]/Ab[i]
			Jx[i] = Kp[i]/Klemmug
			g0[i] = (g[i]/Ap[i])/(Rho_P*rp0[i])
			Re[i] = ((Rho_P*rp0[i]*4*Ap[i])/S[i])/Viscosidad_camara
			gRe[i] = g0[i]*(Re[i]/1000)**(-1/8)
		   
			eta_temp[i] =	1 + 0.023*(g0[i]**0.8-gth**0.8)  -- No funciona sin una variable intermedia
			eta[i] = IF(eta_temp[i]<1) 1 ELSE eta_temp[i]

			rp0[i] = Factor_a*P[i]**Exponente
			
			-- Remanso
			Rho[i] = P[i] /(R_gas*T[i])
			SoundSpeed[i]= sqrt(gamma*R_gas*T[i])
			MACH[i]= U[i]/SoundSpeed[i]
			Tt[i] = T[i]*(1+0.5*(gamma-1)*MACH[i]**2)
			Pt[i] = P[i]*(Tt[i]/T[i])**(gamma/(gamma-1))
		   S[i]'  =	2*Pi*   rp0[i]*eta[i]
			Ap[i]' =  S[i]*  rp0[i]*eta[i]
		END EXPAND_BLOCK
		
		EXPAND_BLOCK (i IN 1,nodos-1)		
		   Ab[i+1] = Ab[i]+0.5*(S[i]+S[i+1])*dx
			
			dg[i+1] = 0.5* Rho_P*(rp0[i]*S[i]+rp0[i+1]*S[i+1])*dx  
			
			g[i+1] = g[i]+ dg[i+1]   
			
			P[i+1] = P[i]-((Rho[i]*U[i]*Ap[i])*(U[i+1]-U[i])+g[i+1]*U[i+1])/(0.5*(Ap[i+1]+Ap[i]))
			
			U[i+1] = ((Rho[i]*U[i]*Ap[i])+ dg[i+1])/(Rho[i+1]*Ap[i+1]) 
			
			T[i+1] = T[i] - ((0.5*(gamma-1))/(gamma*R_gas))  *  (U[i+1]**2 - U[i]**2) 
			
			Coord[i+1] = Coord[i] + dx	
		END EXPAND_BLOCK	
	

END COMPONENT




	
	

