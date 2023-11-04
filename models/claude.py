from dotenv import load_dotenv

import os
import yaml

if __name__ == "__main__":
    load_dotenv()

    with open('config.yml', 'r') as file:
        data = yaml.safe_load(file)

    api_key = os.getenv("CLAUDE_SECRET")