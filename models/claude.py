from dotenv import load_dotenv
from anthropic import Anthropic, HUMAN_PROMPT, AI_PROMPT

import os
import yaml


class Claude(Anthropic):
    def __init__(self):
        self.config = self.config()

        super().__init__(api_key=self.get_api_key())

    def get_api_key(self):
        load_dotenv(dotenv_path=self.get_api_path())
        api_key = os.getenv("CLAUDE_SECRET")
        return api_key

    def get_api_path(self):
        api_path = self.config['paths']['secret']
        return api_path

    @staticmethod
    def config():
        with open('../configs/claude.yml', 'r') as file:
            config = yaml.safe_load(file)
        return config

    def create_completion(self, prompt):
        completion = self.completions.create(
            model="claude-2",
            max_tokens_to_sample=300,
            prompt=prompt
        )
        return completion.completion


if __name__ == "__main__":
    claude = Claude()
    prompt = f"{HUMAN_PROMPT} What is the best way to embed text into an embedding vector? {AI_PROMPT}"
    response = claude.create_completion(prompt)
    print(response)
