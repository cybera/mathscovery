import os
import toml

_this_dir = os.path.dirname(os.path.realpath(__file__))
_settings_file = os.path.join(_this_dir, "../project-settings.toml")

settings = toml.load(_settings_file)
