function () {
  tabPanel ("About", 
            p(style="text-align:justify", ''),
            br(),
            strong('Author'),
            p('Jui-Hua Hsieh',br(),
              a('LinkedIn', href="http://www.linkedin.com/pub/jui-hua-hsieh/12/800/558", target="_blank")
            ),
            br(),
            strong('Assay Descriptions'),
            p(style="text-align:justify", 'Please check descriptions in PubChem', 
              a('PubChem:Tox21', href="https://www.ncbi.nlm.nih.gov/pcassay?term=%22tox21%22", target="_blank")
            ),
            br(),
            strong('Acronyms'),
            p(style="text-align:justify", 
              em('ratio'), ' and ', em('luc'), ' main readout in either beta-lactamase assay or luciferase assay, respectively.', br(),
              em('ch1'), ' indicates the channel 1, the background in bla -  beta-lactamase assay - assay.',  br(), 
              em('ch2'), ' indicates channel 2, the signal channel in bla assay.', br(), 
              em('via'), ' indicates cell viability.', br(),
              em('p-length'), ' and ', em('f-length'), ' indicates partial length (in HEK293 cell) or full length receptor.',
              'Numbers represent different batches.'
              ),
            br()
  )  
}