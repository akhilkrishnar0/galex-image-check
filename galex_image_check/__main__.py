import argparse
import pandas as pd
from .core import CheckConfig,run_check

def main():
    p=argparse.ArgumentParser(description='Check local GALEX intbgsub field coverage')
    p.add_argument('catalog',help='CSV with galaxy,ra_deg,dec_deg columns')
    p.add_argument('--base-dir',default='downloaded_galex')
    p.add_argument('--output-dir',default='galex_image_check_results')
    p.add_argument('--radius-arcmin',type=float,default=2.5)
    p.add_argument('--plot-arcmin',type=float,default=3.0)
    p.add_argument('--min-coverage',type=float,default=.95)
    p.add_argument('--allow-zero',action='store_true')
    args=p.parse_args()
    config=CheckConfig(args.radius_arcmin,args.plot_arcmin,args.min_coverage,not args.allow_zero)
    table=run_check(pd.read_csv(args.catalog),args.base_dir,args.output_dir,config)
    print(table.groupby(['galaxy','band','status']).size().rename('count').to_string())
    print(f'TOTAL: {len(table)} | GOOD: {sum(table.status=="OK")} | FLAGGED: {sum(table.status!="OK")}')

if __name__=='__main__':
    main()
