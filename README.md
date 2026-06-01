# drug-target-repo
## Background

A bioinformatician in our company (who can work with simple APIs and CLIs but not so much with wrangling dataframes/tables)
wants to explore OpenTargets data about drugs and their targets. Specifically,
they want to know for a given new target of interest what clinical or marketed compounds already exist
acting on that target. As there might be many drugs associated with a given target, our users
would like to be pointed to the drugs that are most specific to this particular target (i.e., have associations with as few other targets as possible)
and whose reported adverse effects have the least significant associations to this drug (by [LLR](https://openfda.shinyapps.io/LRTest/_w_c5c2d04d/lrtmethod.pdf), pick a simple, yet reasonable aggregation).
Unfortunately, the information about drug-target associations only lists Ensembl gene identifiers, which means you will need
to map from EnsemblID to Hugo gene symbol (as this is what the biologists will understand).
We recommend the [biothings client](https://github.com/biothings/biothings_client.py) if you are using python.
You might encounter small formatting differences in other identifiers as well, which should be easy to fix by yourself.

## Implementation

You are free to choose your technology stack.

- Provide the user with an API or CLI for easy retrieval of summarized data needed to answer their question described in the Background.
- Impersonate the user and write a small script that uses your API to retrieve the data and suggest a method to pick the best drugs (e.g., a small visualization or scoring)

A potential solution could involve a simple ETL script that fills a database
plus a simple RestAPI that queries that database. Your user script then interacts with the RestAPI.
But if you want to choose a different route and/or feel more comfortable with a different setup,
this is also fine, as long as it gets the job done!
In any case, however, we would like to ask you to apply good software engineer practices while developing your solution.

Ready? Here is your data!

## Data

Data is located in the `./data` folder as gzipped tsvs.

- adverseEffects.tsv.gz contains data about different adverse effects reported by patients for a given drug
- molecules.tsv.gz contains data about the drugs/molecules (and their associations to targets and diseases)

The columns should be rather self-explanatory but if you have any questions, please write us. Asking questions
is encouraged.

## Time scope

We expect the challenge to take around 3 hours of focussed work. If you feel you will exceed this time frame
by a lot, consider resorting to simplifications instead of spending too much time or let alone giving up.
We will focus on your thought process and development practices.

## Solution
Please checkout: [solution.md](https://github.com/Arslan-Siraj/drug-target-repo/blob/main/solution.md)