/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: MChSolido
 AUTHOR: Álvaro González Villarreal
 DESCRIPTION: Seccion de motor cohete de propelente solido
 CREATION DATE: 19/03/2022
-----------------------------------------------------------------------------------------*/
COMPONENT MChSolido
	-- De momento sin puertos	
	DATA
			-- DATOS GENERALES COMBUSTIBLE Y PARAMETROS COMBUSTION
			REAL Rho_P = 1750	UNITS "kg/m3"
			REAL gamma =	1.33	
			REAL Peso_molecular =	27	UNITS "g/mol"
			REAL G_gamma= 0.6726	
			REAL Viscosidad_camara =	1.80E-05	
			BOOLEAN Combustion_erosiva = TRUE	
			REAL Comb_Erosiva_gth =	35	
			REAL Temperatura_de_combustion =3300 UNITS	"K"
			REAL ce_estrella =	1498.7	UNITS "m/s"
			REAL Velocidad_recesion_70 =	7	UNITS "mm/s"
			REAL Exponente =	0.4	
			REAL Factor_a =	1.27955E-05	
			REAL R_gas = 307.93  UNITS "J/kgK"			
			-- DATOS MOTOR
			REAL Area_de_garganta	 = 8.00E-03	UNITS "m2"
			REAL Presion_de_camara= 4.20E+06	UNITS "Pa"
			REAL Tiempo_de_combustion = 85	UNITS "s"
			REAL Klemmug = 280.64	
			REAL Area_de_Combustion =	2.25E+00 UNITS "m2"
			REAL Parametro_Jota_inicial=0.6	
			REAL Num_Mach_port	= 0.350
			REAL Par_func_Y =	0.404
			REAL Area_port = 1.33E-02 UNITS "m2"
			REAL Diametro_port = 1.30E-01 UNITS "m" 
			REAL Perimetro_port = 4.09E-01 UNITS "m" 
			REAL LD_port	= 42.10
			REAL web = 5.95E-01 UNITS "m" 
			REAL Diametro_camara =	1.32E+00 UNITS "m" 
			REAL Gasto_f =	2.24E+01 UNITS "kg/s" 
			REAL Coef_Volumetrico =	0.990
			REAL Esbeltez_LD_camara =	4.15
			REAL L	= 5.48E+00 UNITS "m" 
			-- DATOS SIM
			REAL nodos = 11
			
			
			REAL S[11] = {4.09E-01,4.09E-01,4.09E-01,4.09E-01,4.09E-01,4.09E-01,4.09E-01,4.09E-01,4.09E-01,4.09E-01,4.09E-01} UNITS "m"
			REAL Ap[11]= {1.33E-02,1.33E-02,1.33E-02,1.33E-02,1.33E-02,1.33E-02,1.33E-02,1.33E-02,1.33E-02,1.33E-02,1.33E-02} UNITS "m2"
				
	DECLS
			REAL Coord [11]  UNITS "m"
					
			REAL dx
						
			REAL Ab[11] UNITS "m2"
			--REAL S[11] UNITS "m"
			--REAL Ap[11] UNITS "m2"
			REAL Kp[11]
			REAL Jx[11]
			REAL g0[11]  UNITS "no_units"
			REAL gRe[11] UNITS "no_units"
			REAL Re[11]
			REAL eta[11]
			REAL g[11]
			REAL rp0[11] UNITS "no_units"	--(apc^n)
			REAL dg[11]
			REAL P[11] UNITS "Pa"
			REAL U[11] UNITS "m/s"
			REAL T[11] UNITS "K"
			REAL Rho[11] UNITS "kg/m3"
			REAL SoundSpeed[11]
			REAL MACH[11]
			REAL Tt[11]
			REAL Pt[11]
	
	

		
	CONTINUOUS	
		-- Condiciones de contorno
		dx = L/(nodos-1)
		Coord[1]=0
		Ab[1] = 0
		--dg[1] = 0
		g[1] = 0
		U[1]= 0
		T[1] = Temperatura_de_combustion
		P[1] = 6.33E+06 --( P[1]*(g[10]*ce_estrella))/(Area_de_garganta*Pt[10])
		
		EXPAND_BLOCK (i IN 1,11)
			SoundSpeed[i]= sqrt(gamma*R_gas*T[i])
			MACH[i]= U[i]/SoundSpeed[i]
			Kp[i] = Ap[i]/Ab[i]
			Jx[i] = Kp[i]/Klemmug
			g0[i] = (g[i]/Ap[i])/(Rho_P*rp0[i])
			Re[i] = ((Rho_P*rp0[i]*4*Ap[i])/S[i])/Viscosidad_camara
			gRe[i] = g0[i]*(Re[i]/1000)**(-1/8)
		   eta[i] = 1 + 0.023*(g0[i]**0.8-Comb_Erosiva_gth**0.8)	
			rp0[i] = Factor_a*P[i]**Exponente
			Rho[i] = P[i] /(R_gas*T[i])
			Tt[i] = T[i]*(1+0.5*(gamma-1)*MACH[i]**2)
			Pt[i] = P[i]*(Tt[i]/T[i])**(gamma/(gamma-1))
		END EXPAND_BLOCK
		EXPAND_BLOCK (i IN 1,10)
			-- Falta meter tiempo Leyes temporales de momento solo t=0
			--S[]  = 4.09E-01	-- S(to)+2*PI()*rp0(i)*eta[i]*tiempo
			--Ap[] = 1.33E-02  --Ap(To)+S[i]*eta[i]*rp0[i]*tiempo
			
			Ab[i+1] = Ab[i]+(S[i]+S[i+1])*dx
			dg [i+1] = 0.5* Rho_P*(rp0[i]*S[i]+rp0[i+1]*S[i+1])*dx --Añadir terminos de masa de ignicion?
			g[i+1] = g[i]+ dg[i+1]   
			U[i+1] = ((Rho[i]*U[i]*Ap[i])+dg[i+1] )/(Rho[i+1]*Ap[i+1])
			P[i+1] = P[i] - ((Rho[i]*U[i]*Ap[i])*(U[i+1]-U[i])+g[i+1]*U[i+1])/(0.5*(Ap[i+1]+Ap[i]))
			T[i+1] = T[i] - (0.5*(gamma-1)/(gamma*R_gas))*(U[i+1]**2-U[i]**2)
			
			-- Variacion de rp0?
			
			Coord[i+1] = Coord[i]+ dx	
		END EXPAND_BLOCK	
		
	

END COMPONENT




	
	

