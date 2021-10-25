import fO2calculate as f

my_samples = ['gasdata/SampleA',
			  'gasdata/SampleB',
			  'gasdata/SampleC',
			  'gasdata/SampleD',
			  'gasdata/SampleE',
			  'gasdata/SampleF',
			  'gasdata/SampleG',
			  'gasdata/SampleH',
			  'gasdata/SampleI',
			  'gasdata/SampleJ',
			  'gasdata/SampleK',
			  'gasdata/SampleL',
			  'gasdata/SampleM',
			  'gasdata/SampleN',
			  'gasdata/SampleO']

my_sample_names = ['Sample A',
				   'Sample B',
				   'Sample C',
				   'Sample D',
				   'Sample E',
				   'Sample F',
				   'Sample G',
				   'Sample H',
				   'Sample I',
				   'Sample J',
				   'Sample K',
				   'Sample L',
				   'Sample M',
				   'Sample N',
				   'Sample O']

my_conditions = ['Sil-1',
				 'Sil-1',
				 'Sil-1',
				 'Sil-1',
				 'Sil-2',
				 'Sil-2',
				 'Sil-2',
				 'Sil-2',
				 'Sil-3',
				 'Sil-3',
				 'Sil-3',
				 'Sil-3',
				 'Sil-2',
				 'Sil-3',
				 'Sil-3']

for i, samp in enumerate(my_samples):
	myfluid = f.Fluid(samp, highT=False)
	fO2_C = myfluid.calc_fO2_from_gas_ratios('CO2')
	fO2_H = myfluid.calc_fO2_from_gas_ratios('H2O')
	fCO = myfluid.calc_fCO()
	# make_plots
	f.make_plots(myfluid, 
				 species=['H2', 'CH4', 'H2O', 'CO', 'O2', 'CO2'],
				 fCO=fCO, fO2_from_C_ratios=fO2_C,
				 fO2_from_H_ratios=fO2_H, sample_name=my_sample_names[i],
				 conditions=my_conditions[i])

	# make comparison plots
	f.make_comparison_plots(myfluid, sample_name=my_sample_names[i],
							conditions=my_conditions[i])

