require('torch')
require('cutorch')
require('nn')
require('nnx')
require('cunn')
require('cunnx')
require('hdf5')
require('lfs')
require('io')
require('os')

package.path = package.path .. ";" .. lfs.currentdir() .. '/../../../?/th/init.lua'
lib = require('lib')

-- Private functions definition ------------------------------------------------
local function nilling(module, fields)
   for key, val in pairs(module) do
      if string.match(torch.typename(val) or '', 'Tensor') and not fields[key] then
         module[key] = torch.Tensor():typeAs(val)
      end
   end
end

local function netLighter(network, fields)
   nilling(network, fields)
   if network.modules then
      for _,a in ipairs(network.modules) do
         netLighter(a, fields)
      end
   end
end

local function craftGradNBias(module)
   if module.weight then module.gradWeight = module.weight:clone() end
   if module.bias   then module.gradBias   = module.bias  :clone() end
   if module.__typename == 'nn.SpatialConvolutionCUDA' then
      module.gradWeightPartial = module.weight:clone()
   end
   if module.__typename == 'nn.SpatialConvolutionMM' then
      module.fgradInput = torch.Tensor():typeAs(module.output)
   end
   if type(module.output) == 'table' then
      module.gradInput = {}
      for i, _ in ipairs(module.output) do
          module.gradInput[i] = torch.Tensor():typeAs(module.output[i])
      end
   else
      module.gradInput = torch.Tensor():typeAs(module.output)
   end
end

local function repopulateGradNBias(network)
   craftGradNBias(network)
   if network.modules then
      for _,a in ipairs(network.modules) do
         repopulateGradNBias(a)
      end
   end
end

-- Public functions definition -------------------------------------------------
local function saveLight(fileName, model, fields, format)
   -- Reverse dictionary
   local keepFields = {}
   for _, val in pairs(fields) do
      keepFields[val] = true
   end
   -- Getting rid of unnecessary things and freeing the memory
   netLighter(model, keepFields)
   collectgarbage()
   fileName = format == 'ascii' and fileName .. '.ascii' or fileName
   torch.save(fileName, model, format)
   -- Repopulate the gradWeight through the whole net
   repopulateGradNBias(model)
   -- Return NEW storage for <weight> and <grad>
   return model:getParameters()
end

local function loadLight(fileName, format)
   fileName = format == 'ascii' and fileName .. '.ascii' or fileName
   local model = torch.load(fileName, format)
   -- Repopulate the gradWeight through the whole net
   repopulateGradNBias(model)
   return model--, model:getParameters()
end

function zeroDataSize(data)
  if type(data) == 'table' then
    for i = 1, #data do
      data[i] = zeroDataSize(data[i])
    end
  elseif type(data) == 'userdata' then
    data = torch.Tensor():typeAs(data)
  end
  return data
end

-- Resize the output, gradInput, etc temporary tensors to zero (so that the
-- on disk size is smaller)
function cleanupModel(node)
  if node.output ~= nil then
    node.output = zeroDataSize(node.output)
  end
  if node.gradInput ~= nil then
    node.gradInput = zeroDataSize(node.gradInput)
  end
  if node.finput ~= nil then
    node.finput = zeroDataSize(node.finput)
  end
  -- Recurse on nodes with 'modules'
  if (node.modules ~= nil) then
    if (type(node.modules) == 'table') then
      for i = 1, #node.modules do
        local child = node.modules[i]
        cleanupModel(child)
      end
    end
  end

  collectgarbage()
end

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

if not fileExists(arg[1]) then
  print('Input file does not exist.')
  os.exit()
end

if fileExists(arg[2]) then
  print('Output file already exists. Press ANY key to overwrite.')
  io.read()
end

model = torch.load(arg[1])
print('Read ' .. arg[1] .. '.')
print(model)
cleanupModel(model)
torch.save(arg[2], model)
print('Wrote ' .. arg[2] .. '.')

-- Test
--model = torch.load(arg[2])
--input = torch.Tensor(16, 2, 24, 54, 24):fill(0)
--input = input:cuda()
--model:forward(input)