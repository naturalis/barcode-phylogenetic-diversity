##Tried something like this. Says: The AWS Access Key Id you provided does not exist in our records.
import boto3
import os
from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv(usecwd=True))

aws_access_key = os.environ.get('AWS_ACCESS_KEY_ID')
aws_secret_key = os.environ.get('AWS_SECRET_ACCESS_KEY')

s3 = boto3.client('s3',\
                  aws_access_key_id=aws_access_key,\
                  aws_secret_access_key=aws_secret_key,\
                  )
list_buckets = s3.list_buckets()