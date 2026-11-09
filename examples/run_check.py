from pathlib import Path
import pandas as pd
from galex_image_check import run_check,build_exclusions
catalog=pd.read_csv(Path(__file__).with_name('catalog.csv'))
df=run_check(catalog)
print(df.groupby(['galaxy','band','status']).size())
print(build_exclusions(df))
