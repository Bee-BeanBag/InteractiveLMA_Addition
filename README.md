# InteractiveLMA_Addition
An addition of functions to enhance pyxlma by deeplycloudy
Currently tested for windows. Linux and Mac ports in development


|Installation guide|
```
git clone https://github.com/Bee-BeanBag/InteractiveLMA_Addition.git
cd InteractiveLMA_Addition
conda env create -f environment.yml

```

From here you are able to make use of all the additions! Examples of the functions can be found in test_code.py with the source for each function in Thesis_Code_Functs.py

|Uses|
- Graph speed of lightning leaders relative to first LMA source
- Plot LMA sources over radar plan view and cross section
- Create gridded average base reflectivity plot for a list of radar scans
- Create base radar frequency plot for list of radar scans
- Create flash extent density for all the sources plotted on the current LMA plot
    - Combine last three into one plot 
