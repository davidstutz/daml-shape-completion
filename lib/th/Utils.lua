-- Some utilities.

-- https://github.com/harningt/luajson
require('json')
-- https://github.com/deepmind/torch-hdf5
require('hdf5')
-- http://keplerproject.github.io/luafilesystem
require('lfs')

--- @module utils
local utils = {}

--- Recursively prints a table and all its subtables.
-- @see https://coronalabs.com/blog/2014/09/02/tutorial-printing-table-contents/
-- @param t table to print
function utils.printTable(t)

  -- A cache for all printed tables.
  local printCache = {}

  local function subPrintTable(t, indent)
    if (printCache[tostring(t)]) then
      print(indent .. '*' .. tostring(t))
    else
      printCache[tostring(t)]=true
      if (type(t) == 'table') then
        for pos,val in pairs(t) do
          if (type(val) == 'table') then
            print(indent .. '[' .. pos .. '] => ' .. tostring(t) .. ' {')
            subPrintTable(val, indent..string.rep(' ', string.len(pos) + 8))
            print(indent .. string.rep(' ', string.len(pos) + 6) .. '}')
          elseif (type(val) == 'string') then
            print(indent .. '[' .. pos .. '] => "' .. val .. '"')
          else
            print(indent .. '[' .. pos .. '] => ' .. tostring(val))
          end
        end
      else
        print(indent .. tostring(t))
      end
    end
  end

  if (type(t) == 'table') then
    print(tostring(t) .. ' {')
    subPrintTable(t, '  ')
    print('}')
  else
    subPrintTable(t, '  ')
  end
end

--- Checks if a file exists.
-- @see http://stackoverflow.com/questions/4990990/lua-check-if-a-file-exists
-- @param filePath path to file
-- @return true if file exists
function utils.fileExists(filePath)
  local f = io.open(filePath, 'r')
  if f ~= nill then
    io.close(f)
    return true
  else
    return false
  end
end

--- Checks if a directory exists using the lfs package.
-- @param dirPath path to directory
-- @return true if directory exists
function utils.directoryExists(dirPath)
  local attr = lfs.attributes(dirPath)
  if attr then
    if attr['mode'] == 'directory' then
      return true
    end
  end

  return false
end

--- Reverse a list.
-- @see http://lua-users.org/wiki/ListOperations
-- @param list list ot reverse
-- @return reversed list
function utils.reverseList(list)
  local rList = {}
  for i = table.getn(list), 1, -1 do
    table.insert(rList, list[i])
  end
  return rList
end

--- Recursively create the given directory; not throrougly tested, might be sensitive to non-linux
-- file paths.
-- @param dirPath path to directory
function utils.makeDirectory(dirPath)
  local function findDirectories(subPath, dirCache, i)
    local lastChar = dirPath:sub(subPath:len(), subPath:len())
    if lastChar == '/' then
      subPath = subPath:sub(1, -2)
    end

    if subPath:len() > 0 then
      if not utils.directoryExists(subPath) then
        dirCache[i] = subPath
        -- http://stackoverflow.com/questions/5243179/what-is-the-neatest-way-to-split-out-a-path-name-into-its-components-in-lua
        local subSubPath, subDir, ext = string.match(subPath, "(.-)([^\\/]-%.?([^%.\\/]*))$")
        findDirectories(subSubPath, dirCache, i + 1)
      end
    end
  end

  local dirCache = {}
  findDirectories(dirPath, dirCache, 1)
  local rDirCache = utils.reverseList(dirCache)

  for i = 1, #rDirCache do
    lfs.mkdir(rDirCache[i])
  end
end

--- Took me 20 minutes to figure out that LUA/Torch are so f***ing stupid that this
-- is not possible without iterating!
-- @param storage storage to compute product of
-- @return product of all dimensions
function utils.storageProd(storage)
  if #storage == 0 then
    return 0
  end

  local prod = 1
  for i = 1, #storage do
    prod = prod * storage[i]
  end
  return prod
end

--- Compute the sum of storage elements.
-- @param storage storage to compute product of
-- @return product of all dimensions
function utils.storageSum(storage)
  local sum = 0
  for i = 1, #storage do
    sum = sum + storage[i]
  end
  return sum
end

--- Write a table as JSON to a file.
-- @param file file to write
-- @param t table to write
function utils.writeJSON(file, t)
  local f = assert(io.open(file, 'w'))
  f:write(json.encode(t))
  f:close()
end

--- Read a JSON file into a table.
-- @param file file to read
-- @return JSON string
function utils.readJSON(file)
  local f = assert(io.open(file, 'r'))
  local tJSON = f:read('*all')
  f:close()
  return json.decode(tJSON)
end

--- Writes a single torch tensor to HDF5.
-- @param file file to write to
-- @param tensor tensor to write
-- @param key optional key, i.e. tensor is accessible as "/key"
function utils.writeHDF5(file, tensor, key)
  local key = key or 'tensor'
  local h5 = hdf5.open(file, 'w')
  h5:write('/' .. key, tensor)
  h5:close()
end

--- Reads a single torch tensor from HDF5.
-- @param file file to read
-- @param key key to read from, i.e. read "/key"
-- @return tensor
function utils.readHDF5(file, key)
  local key = key or 'tensor'
  local h5 = hdf5.open(file, 'r')
  tensor = h5:read('/' .. key):all()
  h5:close()
  return tensor
end

--- Copies the weights of the given layers between two models; assumes the layers to have .weight and .bias defined.
-- @param modelFrom mode to copy weights from
-- @param modelTo model to copy weights to
-- @param layersFrom layer indices in model_from
-- @param layersTo layer indices in model_to
function utils.copyWeights(modelFrom, modelTo, layersFrom, layersTo)
  assert(#layersFrom == #layersTo)

  for i = 1, #layersFrom do
    --if modelTo.modules[layersTo[i]].weight ~= nil or modelTo.modules[layersTo[i]].bias ~= nil then
      assert(modelFrom.modules[layersFrom[i]].__typename == modelTo.modules[layersTo[i]].__typename,
          'layer from ' .. layersFrom[i] .. ' and layer to ' .. layersTo[i] .. ' are not of the same type!')

      -- Allows to provide all layers, also these without parameters.
      if modelTo.modules[layersTo[i]].weight ~= nil then
        modelTo.modules[layersTo[i]].weight = modelFrom.modules[layersFrom[i]].weight:clone()
        modelTo.modules[layersTo[i]].gradWeight:resize(#modelFrom.modules[layersFrom[i]].gradWeight)
      end
      if modelTo.modules[layersTo[i]].bias ~= nil then
        modelTo.modules[layersTo[i]].bias = modelFrom.modules[layersFrom[i]].bias:clone()
        modelTo.modules[layersTo[i]].gradBias:resize(#modelFrom.modules[layersFrom[i]].gradBias)
      end
    --end
  end
end

--- Copies the weights to a subnetwork. The subnetwork is expected to have the same
-- structure and optionally start at the provided layer index.
-- @param modelFrom model to copy weights from
-- @param modelTo model to copy weights to; expected to be a subnetwork starting at startLayer
-- @param fromStart start layer in modelFrom
-- @param toStart start layer in modelTo
-- @param numLayers number of layers
function utils.copyWeightsSubNetwork(modelFrom, modelTo, fromStart, toStart, numLayers)
  fromStart = fromStart or 1
  toStart = toStart or 1
  numLayers = numLayers or math.min(#modelFrom.modules - fromStart + 1, #modelTo.modules - toStart + 1)

  local layersFrom = {}
  local layersTo = {}
  for i = 1, numLayers do
    layersFrom[i] = (fromStart - 1) + i
    layersTo[i] = (toStart - 1) + i
  end

  --print(modelFrom)
  --print(modelTo)
  --print(layersFrom)
  --print(layersTo)

  utils.copyWeights(modelFrom, modelTo, layersFrom, layersTo)
end

--- Sets all layers with parameters (weights or biases) to be fixed, i.e. overwrites
-- the paramters function to return nothing and the accGradParameters function to
-- to nothing. Should be applied before getParameters is called!
-- @param model model to fix the given layers
-- @param layers indices of layers to fix.
function utils.fixLayers(model, layers)
  for i = 1, #layers do
    if model.modules[layers[i]].weight ~= nil or model.modules[layers[i]].bias ~= nil then

      -- Set gradients to nil for clarity.
      if model.modules[layers[i]].weight ~= nil then
        model.modules[layers[i]].gradWeight = nil
      end
      if model.modules[layers[i]].bias ~= nil then
        model.modules[layers[i]].gradBias = nil
      end

      -- Has no trainable parameters.
      model.modules[layers[i]].parameters = function() end
      -- Does not compute gradients w.r.t. parameters.
      model.modules[layers[i]].accGradParameters = function(input, gradOutput, scale) assert(model.modules[layers[i]].gradWeight == nil) end
      -- Note that updateGradInput is not touched!
    end
  end
end

--- Sets all layers with parameters (weights and biases) to be fixed starting with the given
-- start layers.
-- @param model model to fix layers
-- @param startLayer starting layer
function utils.fixLayersAfter(model, startLayer)
  local j = 1
  local layers = {}

  for i = startLayer, #model.modules do
    layers[j] = i
    j = j + 1
  end

  utils.fixLayers(model, layers)
end

--- Find all layers of the given type.
-- @param model model to look in
-- @param type type name of the layers to look for
-- @return layers in order
function utils.findLayers(model, type)
  local j = 1
  local layers = {}

  for i = 1, #model.modules do
    if model.modules[i].__typename == type then
      layers[j] = model.modules[i]
      j = j + 1
    elseif model.modules[i].mdoules ~= nil then
      local subLayers = utils.findLayers(model.modules[i], type)
      for k = 1, #subLayers do
        layers[j] = subLayers[k]
        j = j + 1
      end
    end
  end

  return layers
end

--- Finds the first layer of the given type.
-- @param model model to look in
-- @param type type name of the layers to look for
-- @return layer
function utils.findLayerFirst(model, type)
  local layers = utils.findLayers(model, type)
  assert(#layers > 0)
  return layers[1]
end

--- Split text into a list consisting of the strings in text,
-- separated by strings matching delimiter (which may be a pattern).
-- @see http://lua-users.org/wiki/SplitJoin
-- @param delimited delimited to split string by
-- @param text text to split
-- @return table of strings
function utils.splitString(delimiter, text)
  local strfind = string.find
  local strsub = string.sub
  local tinsert = table.insert

  local list = {}
  local pos = 1

  if strfind('', delimiter, 1) then -- this would result in endless loops
    assert(false, 'delimiter matches empty string!')
  end

  while 1 do
    local first, last = strfind(text, delimiter, pos)
    if first then -- found?
      tinsert(list, strsub(text, pos, first-1))
         pos = last+1
    else
      tinsert(list, strsub(text, pos))
      break
    end
  end

  return list
end

--- Append the tensor tensor to the tensor acc which may initially be nil;
-- the given tensor will be converted to float.
-- @param acc accumulated tensor
-- @param tensor tensor to accumulate
function utils.appendTensor(acc, tensor)
  if acc == nil then
    acc = tensor:float()
  else
    acc = torch.cat(acc, tensor:float(), 1)
  end

  return acc
end

lib.utils = utils