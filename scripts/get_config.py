import configparser

def _cast(value):
    value = value.split('#')[0].strip()

    if '%' in value:
        return value

    if value == 'None':
        return None
    try:
        if value.find('.') == -1 and value.find('e') == -1:
            return int(value)
        else:
            return float(value)
    except ValueError:
        # Check if it can be cast as a boolean.
        if value in 'True False'.split():
            return eval(value)
        # Return as the original string.
        return value

def get_config(config_file):
    cp = configparser.ConfigParser()
    cp.optionxform = str
    cp.read(config_file)
    return {k: _cast(v) for k, v in cp.items('DEFAULT')}

if __name__ == '__main__':
    config = get_config('instcat_validation_config.ini')
    print(config)
