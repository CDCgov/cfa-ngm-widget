# Next generation matrix widget

⚠️ This is a work in progress

## Overview

This repo contains code for investigating the potential efficacy of vaccination allocation for a disease interactively via a [streamlit](https://streamlit.io/) app using a next generation matrix approach.

## Getting started

- Enable poetry with `poetry install`
- To run the app: `make`, which calls `streamlit run ngm/widget.py`

### Using containers

The build and run process can also be executed using the `build_container` and `run_container` targets in the included [Makefile](Makefile).

## Model Description

This repo contains code to apply the next-generation method of Diekman et al. (1990) to calculate R0 for an SIR model with 3 risk groups and flexible inputs for varying vaccine allocation to each group.

### Next generation matrix calculation

The next-generation matrix approach (NGM) is described in `docs/ngm.md`. The documentation is best viewed off of GitHub, either by opening in VSCode and using the built in [markdown preview](https://code.visualstudio.com/Docs/languages/markdown#_markdown-preview) or by building with `mkdocs` using `mkdocs serve` (this requires installing `mkdocs`).

### Model assumptions

Vaccination is assumed to be all or nothing -- each individual's immunity is determined by a coin flip with probability of being immune equal to the vaccine efficacy.

### Widget

The widget is designed to let users modify vaccine allocation choices, measures of between- and within-group spread, population composition, and see how this effects $R_e$ and other key quantities.

Specific inputs to and outputs from the widget are documented therein.

### References

Diekmann O, Heesterbeek JA, Metz JA. On the definition and the computation of the basic reproduction ratio R0 in models for infectious diseases in heterogeneous populations. J Math Biol. 1990;28(4):365-82. doi: 10.1007/BF00178324. PMID: 2117040.

Diekmann O, Heesterbeek JA, Roberts MG. The construction of next-generation matrices for compartmental epidemic models. J R Soc Interface. 2010 Jun 6;7(47):873-85. doi: 10.1098/rsif.2009.0386. Epub 2009 Nov 5. PMID: 19892718; PMCID: PMC2871801.

van den Driessche P, Watmough J. Reproduction numbers and sub-threshold endemic equilibria for compartmental models of disease transmission. Math Biosci. 2002 Nov-Dec;180:29-48. doi: 10.1016/s0025-5564(02)00108-6. PMID: 12387915.

## Authors

- Paige Miller <yub1@cdc.gov>
- Scott Olesen <ulp7@cdc.gov>
- Andy Magee <rzg0@cdc.gov>

## General Disclaimer

This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm). GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise.

## Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/). All contributions to this repository will be released under the CC0 dedication. By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

## License Standard Notice

This repository is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the [Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md) and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md). For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo) and submitting a pull request. (If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git).) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the [Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice

This repository is not a source of government records but is a copy to increase collaboration and collaborative potential. All government records will be published through the [CDC web site](http://www.cdc.gov).
