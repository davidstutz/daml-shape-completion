require('hdf5')
require('lfs')
require('io')
require('os')

--- Checks if a file exists.
-- @see http://stackoverflow.com/questions/4990990/lua-check-if-a-file-exists
-- @param filePath path to file
-- @return true if file exists
function fileExists(filePath)
  local f = io.open(filePath, 'r')
  if f ~= nill then
    io.close(f)
    return true
  else
    return false
  end
end

--- Writes a single torch tensor to HDF5.
-- @param file file to write to
-- @param tensor tensor to write
-- @param key optional key, i.e. tensor is accessible as "/key"
function writeHDF5(file, tensor, key)
  local key = key or 'tensor'
  local h5 = hdf5.open(file, 'w')
  h5:write('/' .. key, tensor)
  h5:close()
end

--- Reads a single torch tensor from HDF5.
-- @param file file to read
-- @param key key to read from, i.e. read "/key"
-- @return tensor
function readHDF5(file, key)
  local key = key or 'tensor'
  local h5 = hdf5.open(file, 'r')
  tensor = h5:read('/' .. key):all()
  h5:close()
  return tensor
end

if not fileExists(arg[1]) then
  print('Input file does not exist.')
  os.exit()
end

tensor = readHDF5(arg[1])
print('Read ' .. arg[1] .. '.')
print('Shape:')
print(#tensor)