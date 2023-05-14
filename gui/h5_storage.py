import h5py
import numpy as np
#import tracking

dt = h5py.special_dtype(vlen=bytes)

def stringDataset(group, name, data, system=None):
    dset = group.create_dataset(name, (1,), dtype=dt, data=data)
    if system:
        addSystemAttribute(dset, system)
    return dset

def addStringAttribute(dset_or_group, name, data):
    #return dset_or_group.attrs.create(name, np.string_(data)) # , (1,), dtype=dt)
    dset_or_group.attrs[name] = bytes(data, 'utf-8')

def addSystemAttribute(dset_or_group, data):
    return addStringAttribute(dset_or_group, 'system', data)

def add_dataset(group, name, data, system=None, dtype=None):
    if type(data) is str:
        return stringDataset(group, name, data, system)
    else:
        if dtype:
            dset = group.create_dataset(name, data=data, dtype=dtype)
        else:
            dset = group.create_dataset(name, data=data)
        if system:
            addSystemAttribute(dset, system)
        return dset

def saveH5Recursive(h5_filename, data_dict):
    # import here to prevent circular import

    def recurse_save(group, dict_or_data, dict_or_data_name, new_group=None):
        if dict_or_data is None:
            dict_or_data = 'None'
        if group is None:
            print("'recurse_save' has been called with None")
            raise ValueError

        if hasattr(dict_or_data, 'to_dict_custom'):
            dict_or_data = dict_or_data.to_dict_custom()

        if type(dict_or_data) is tuple:
            dict_or_data = {'tuple_%i': x for i, x in enumerate(dict_or_data)}

        if type(dict_or_data) is dict:

            if type(dict_or_data_name) is int:
                inner_key = str(dict_or_data_name)
            else:
                inner_key = dict_or_data_name
            try:
                new_group = group.create_group(inner_key)
            except Exception as e:
                print(e)
                print(dict_or_data_name, 'error')
                #raise
            if new_group is None:
                raise ValueError
            for key, val in dict_or_data.items():
                try:
                    recurse_save(new_group, val, key)
                except ValueError as e:
                    print(e)
                    print('I called recurse_save with None')
                    #import pdb; pdb.set_trace()
        else:
            mydata = dict_or_data
            inner_key = dict_or_data_name
            if type(mydata) is str:
                add_dataset(group, inner_key, mydata.encode('utf-8'), 'unknown')
            #elif type(mydata) is tuple:
            #    mydata2 = np.array(mydata)
            #    add_dataset(group, inner_key, mydata2, 'unknown')
            elif (type(mydata) is list and type(mydata[0]) is str) or (hasattr(mydata, 'dtype') and mydata.dtype.type is np.str_):
                # For list of strings, we need this procedure
                try:
                    if hasattr(mydata, 'dtype') and mydata.dtype.type is np.str and len(mydata.shape) == 2:
                        mydata = mydata.flatten()
                    if len(mydata.shape) == 2:
                        new_list = [[n.encode('ascii') for n in arr] for arr in mydata]
                        max_str_size = max(max(len(n) for n in arr) for arr in mydata)
                    elif len(mydata.shape) == 1:
                        new_list = [n.encode('ascii') for n in mydata]
                        max_str_size = max(len(n) for n in mydata)
                    elif len(mydata.shape) == 0:
                        new_list = [mydata.encode('ascii')]
                        max_str_size = len(new_list[0])
                    #print('Max len %i' % max_str_size)
                    dset = group.create_dataset(inner_key, mydata.shape, 'S%i' % max_str_size, new_list)
                    #print(np.array(dset))
                    dset.attrs.create('system', 'unknown', (1,), dtype=dt)

                except:
                    print('Error for key', inner_key)
                    print(type(mydata))
                    if type(mydata) is list:
                        print('type(mydata[0])')
                        print(type(mydata[0]))
                    print('mydata')
                    print(mydata)

            elif hasattr(mydata, 'dtype') and mydata.dtype == np.dtype('O'):
                if mydata.shape == ():
                    add_dataset(group, inner_key, mydata, 'unknown')
                elif len(mydata.shape) == 1:
                    try:
                        add_dataset(group, inner_key, mydata, 'unknown')
                    except Exception as e:
                        print(e)
                        print('Error for key', inner_key)
                        print(group, inner_key)
                else:
                    for i in range(mydata.shape[0]):
                        for j in range(mydata.shape[1]):
                            try:
                                add_dataset(group, inner_key+'_%i_%i' % (i,j), mydata[i,j], 'unknown')
                            except:
                                print('Error for key', inner_key)
                                print(group, inner_key, i, j)
            elif type(mydata) is list:
                try:
                    add_dataset(group, inner_key, mydata, 'unknown')
                except Exception:
                    mydata2 = {'list_entry_%i' % k: mydata[k] for k in range(len(mydata))}
                    recurse_save(group, mydata2, inner_key)
            else:
                try:
                    add_dataset(group, inner_key, mydata, 'unknown')
                except Exception as e:
                    print('Error', e)
                    print(inner_key, type(mydata))

    with h5py.File(h5_filename, 'w') as dataH5:
        for main_key, subdict in data_dict.items():
            recurse_save(dataH5, subdict, main_key, None)
        #recurse_save(dataH5, data_dict, 'none', new_group=dataH5)

def loadH5Recursive(h5_file):
    def recurse_load(group_or_val, key, saved_dict_curr):
        type_ = type(group_or_val)

        try:
            group_or_val[()]
        except:
            hasval = False
        else:
            hasval = True

        if type_ is h5py._hl.files.File:
            for new_key, new_group_or_val in group_or_val.items():
                recurse_load(new_group_or_val, new_key, saved_dict_curr)
        elif type_ is h5py._hl.group.Group:
            saved_dict_curr[key] = new_dict = {}
            for new_key, new_group_or_val in group_or_val.items():
                recurse_load(new_group_or_val, new_key, new_dict)
        elif type_ == np.dtype('O') and hasval and type(group_or_val[()]) is bytes:
            saved_dict_curr[key] = group_or_val[()].decode()
        elif type_ == h5py._hl.dataset.Dataset:
            dtype = group_or_val.dtype
            #if not hasattr(group_or_val, 'value'):
            #    print('Could not store key %s with type %s in dict' % (key, dtype))
            #    return
            if dtype in (np.dtype('int64'), np.dtype('int32'), np.dtype('int16'), np.dtype('int8'), np.dtype('uint32'), np.dtype('uint16'), np.dtype('uint8'), np.dtype('uint64')):
                saved_dict_curr[key] = np.array(group_or_val[()], dtype).squeeze()
                if saved_dict_curr[key].size == 1:
                    saved_dict_curr[key] = int(saved_dict_curr[key])
            elif dtype == np.dtype('bool'):
                try:
                    saved_dict_curr[key] = bool(group_or_val[()])
                except Exception as e:
                    print(e)
                    print('Could not store key %s with type %s in dict (1)' % (key, dtype))
            elif dtype in (np.dtype('float64'), np.dtype('float32')):
                saved_dict_curr[key] = np.array(group_or_val[()]).squeeze()
                if saved_dict_curr[key].size == 1:
                    saved_dict_curr[key] = float(saved_dict_curr[key])

            elif dtype.str.startswith('|S'):
                if group_or_val[()].shape == (1,1):
                    saved_dict_curr[key] = group_or_val[()][0,0].decode()
                elif group_or_val[()].shape == (1,):
                    saved_dict_curr[key] = group_or_val[()][0].decode()

                elif group_or_val[()].shape == ():
                    saved_dict_curr[key] = group_or_val[()].decode()
                else:
                    saved_dict_curr[key] = [x.decode() for x in group_or_val[()].squeeze()]
            elif dtype.str == '|O':
                saved_dict_curr[key] = group_or_val[()]
            elif type(group_or_val[()]) is str:
                saved_dict_curr[key] = group_or_val[()]
            else:
                print('Could not store key %s with type %s in dict (2)' % (key, dtype))
        else:
            print('Could not store key %s with type %s in dict (3)' % (key, type_))

    saved_dict = {}
    with h5py.File(h5_file.strip(), 'r') as f:
        if 'none' in f:
            recurse_load(f['none'], 'key', saved_dict)
            saved_dict = saved_dict['key']
        else:
            recurse_load(f, 'key', saved_dict)
    return saved_dict

def list_dict_to_list(dd):
    return [dd['list_entry_%i' % i] for i in range(len(dd))]

