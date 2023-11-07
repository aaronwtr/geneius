from dotenv import load_dotenv
from anthropic import Anthropic

import os
import yaml


class Claude(Anthropic):
    def __init__(self, api_key=None):
        self.api_key = api_key

        super().__init__(api_key=self.api_key)

    def create_completion(self, prompt):
        completion = self.completions.create(
            model="claude-2",
            max_tokens_to_sample=600,
            prompt=prompt,
        )
        return completion.completion

    # def get_api_key(self):
    #     load_dotenv(dotenv_path=self.get_api_path())
    #     api_key = os.getenv("CLAUDE_SECRET")
    #     return api_key

    # def get_api_path(self):
    #     api_path = self.config['paths']['secret']
    #     return api_path

    # @staticmethod
    # def config():
    #     with open('configs/claude.yml', 'r') as file:
    #         config = yaml.safe_load(file)
    #     return config
