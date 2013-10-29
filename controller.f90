module Controller
	use trcom

	implicit none

contains

	function computeOutput(omega) result (y)
		real(kind=8), dimension(60), intent(in) :: omega
		real(kind=8), dimension(6) :: y

		y = omega(6:36:6)

	end function computeOutput

	function computeControl(y) result (u)
		real(kind=8), dimension(6), intent(in) :: y
		real(kind=8) :: u
		real(kind=8), dimension(10), save :: xhat = 0
		real(kind=8), dimension(10,10), parameter :: A0 = reshape((/ &
			-8.132180463624117e+04,	5.446682747728348e+04,	1.099458004444427e+05,	-5.519082755715320e+04,	-9.095577501217865e+04,	3.000078953876711e+04,	6.574249098167320e+04,	-1.433961634991070e+04,	-3.800540522829380e+04,	4.381898569011401e+03,	&
			-1.959035875154089e+04,	1.017510472793072e+04,	2.925426464524526e+04,	-1.331474153985595e+04,	-2.542523231826028e+04,	9.325597299735291e+03,	1.605886616772251e+04,	-3.430519304758880e+03,	-1.017990272284520e+04,	2.083724239165059e+03,	&
			8.788740058612571e+03,	-5.519686398446035e+03,	-1.604639199041718e+04,	7.584483491843386e+03,	1.089194763075107e+04,	-5.079616550695963e+03,	-6.326310556899432e+03,	7.830172395365228e+02,	5.120194878844694e+03,	-1.340304699571711e+03,	&
			4.153782997488375e+01,	-5.537957405181292e+02,	8.501777536586457e+02,	-3.590650050874601e+03,	2.537207496506766e+03,	-6.251037235910446e+02,	-1.581362831215129e+03,	1.670970075786345e+03,	-9.873112224828012e+02,	7.817786916075746e+02,	&
			-6.718769455851261e+03,	4.273844085912402e+03,	8.787187482686895e+03,	-2.330356456345257e+03,	-1.293760319375286e+04,	5.845055445569410e+03,	4.711757469003298e+03,	-2.432505931819550e+03,	-1.630667290615545e+03,	-7.820622005597440e+02,	&
			-9.824304076836357e+02,	6.588112975158409e+02,	1.262529364265479e+03,	-1.994644112500017e+03,	2.343883995573714e+03,	-5.905031250346530e+03,	5.146814664651613e+03,	-1.850559478091590e+03,	-1.329712217595457e+03,	1.744415597871003e+03,	&
			2.069962037245166e+03,	-1.398772472907029e+03,	-3.083479833667526e+03,	1.705279195447282e+03,	5.584787441253264e+02,	3.914074544745080e+03,	-9.849723092834502e+03,	6.299757532140172e+03,	-1.526945912979157e+03,	-6.276090088303583e+02,	&
			2.698411248961411e+02,	-3.172624598415323e+02,	-1.770068587990436e+02,	7.723393478819405e+01,	8.316922435526801e+02,	-2.907245747979814e+03,	6.207486724663961e+03,	-1.015106352502444e+04,	7.794451478843444e+03,	-3.665693204688591e+03,	&
			-5.169203680461927e+02,	4.551664895026592e+02,	5.237356989538224e+02,	-1.208487017459576e+02,	-7.769732295091851e+02,	1.242325962459435e+03,	-3.405337389845820e+03,	8.181707225856318e+03,	-1.292896796005763e+04,	9.706295009527805e+03,	&
			7.128221710287156e+01,	-3.708274799299289e+01,	1.102646832382794e+02,	-1.916880213137538e+02,	4.205724628571820e+02,	-3.558287687804864e+02,	1.566817775177927e+03,	-4.973978702368244e+03,	1.039709832186420e+04,	-1.545999307065060e+04	/), shape(A0))
		real(kind=8), dimension(10), parameter :: K = &
			(/ -2.118333808093285e+04,	-5.378691685415339e+03,	2.508578857248183e+03,	7.539212883124128e+01,	-1.836417995747986e+03,	-2.466139931441285e+02,	5.842479508015531e+02,	6.009581584141415e+01,	-1.368372422283719e+02,	1.317528010700335e+01	/)
		real(kind=8), dimension(10,6), parameter :: L = reshape((/ &
			3.114232597523670e+02,	4.236013804382707e+02,	4.106891369918106e+02,	3.619977646561535e+02,	3.131775717687609e+02,	2.153503288955189e+02,	1.296681051962553e+02,	6.791065450790806e+01,	2.646582374441781e+01,	9.133167883847881e-01,	&
			3.626289867049400e+02,	4.156924680464749e+02,	2.173847163326302e+02,	-1.482829459710901e+01,	-1.495416342369639e+02,	-2.246676753134726e+02,	-2.050570253205452e+02,	-1.236399927574198e+02,	-4.146587168682203e+01,	-5.799394017458576e+00,	&
			3.921892533668442e+02,	3.106887968550379e+02,	-9.089446667398003e+01,	-3.352642213230490e+02,	-2.338854233867952e+02,	-3.927793898247338e+01,	9.516587700941395e+01,	1.113190414355605e+02,	5.902768307851707e+01,	7.348827948482684e+00,	&
			3.625486307442808e+02,	1.301076296019734e+02,	-3.048971966862257e+02,	-2.934840308508256e+02,	7.219124146986388e+01,	1.993741789247129e+02,	5.645043233104004e+01,	-8.023844410418420e+01,	-7.052915032388280e+01,	-1.877025334537195e+01,	&
			2.316724440119790e+02,	-8.229726005190942e+00,	-2.081146831575195e+02,	-5.210569450009301e+01,	1.496592057016776e+02,	4.275932500282004e+01,	-7.814643503814457e+01,	-3.321326134232278e+01,	4.346948571469713e+01,	2.630450296456112e+01,	&
			9.880771772153057e+01,	-4.066403108266667e+01,	-4.770882814507560e+01,	2.192076420808040e+01,	2.631430102988362e+01,	-2.475446569035405e+01,	-2.034128279835845e+00,	1.280420343875539e+01,	-6.049443988934261e+00,	-1.292143202846865e+01,	&
			3.114232597523670e+02,	4.236013804382707e+02,	4.106891369918106e+02,	3.619977646561535e+02,	3.131775717687609e+02,	2.153503288955189e+02,	1.296681051962553e+02,	6.791065450790806e+01,	2.646582374441781e+01,	9.133167883847881e-01,	&
			3.626289867049400e+02,	4.156924680464749e+02,	2.173847163326302e+02,	-1.482829459710901e+01,	-1.495416342369639e+02,	-2.246676753134726e+02,	-2.050570253205452e+02,	-1.236399927574198e+02,	-4.146587168682203e+01,	-5.799394017458576e+00,	&
			3.921892533668442e+02,	3.106887968550379e+02,	-9.089446667398003e+01,	-3.352642213230490e+02,	-2.338854233867952e+02,	-3.927793898247338e+01,	9.516587700941395e+01,	1.113190414355605e+02,	5.902768307851707e+01,	7.348827948482684e+00,	&
			3.625486307442808e+02,	1.301076296019734e+02,	-3.048971966862257e+02,	-2.934840308508256e+02,	7.219124146986388e+01,	1.993741789247129e+02,	5.645043233104004e+01,	-8.023844410418420e+01,	-7.052915032388280e+01,	-1.877025334537195e+01	/), shape(L))

		u = -dot_product(K, xhat)
		xhat = xhat + DT0 * (matmul(A0, xhat) + matmul(L, y))

	end function computeControl

	function applyControl(omega, I2) result (torque)
		real(kind=8), dimension(60), intent(in) :: omega
		real(kind=8), intent(in) :: I2
		real(kind=8), dimension(60) :: torque
		real(kind=8), dimension(60), parameter :: g = &
			(/	1.160246547685049e-29,	3.573811599687897e-29,	1.070653826457304e-28,	3.119629144588730e-28,	8.840833768508971e-28,	2.436796977020647e-27,	6.532549211727797e-27,	1.703255782257964e-26,	4.319327021366347e-26,	1.065339810705977e-25,	2.555607234572821e-25,	5.962645963829326e-25,	1.353068633947035e-24,	2.986311043316503e-24,	6.410447582872194e-24,	1.338376090695495e-23,	2.717703225252607e-23,	5.367404605817520e-23,	1.031008733710342e-22,	1.926174588504057e-22,	3.499985058938762e-22,	6.185480006466095e-22,	1.063200534283339e-21,	1.777439763949740e-21,	2.890083774087320e-21,	4.570476977543245e-21,	7.029907284669205e-21,	1.051657237927423e-20,	1.530149668189130e-20,	2.165364531785806e-20,	2.980329807453247e-20,	3.989641153347832e-20,	5.194440406571806e-20,	6.577801527174692e-20,	8.101367988505834e-20,	9.704492869132990e-20,	1.130637319451491e-19,	1.281179776787541e-19,	1.411994875812397e-19,	1.513535250488271e-19,	1.577931240486856e-19,	1.599999999999989e-19,	1.577931616694666e-19,	1.513535370773135e-19,	1.411995044135353e-19,	1.281179471330396e-19,	1.130636420901609e-19,	9.704489398536708e-20,	8.101364608351309e-20,	6.577798390636132e-20,	5.194437620057459e-20,	3.989634811990544e-20,	2.980327853397664e-20,	2.165365564312256e-20,	1.530150458624351e-20,	1.051657822974931e-20,	7.029911474821704e-21,	4.570490780259752e-21,	2.890085726395145e-21,	1.777439446118769e-21	/)
		
		torque = g * I2 * omega

	end function applyControl

end module Controller
