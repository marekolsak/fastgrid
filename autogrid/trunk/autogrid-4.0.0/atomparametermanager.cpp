#include "atomparametermanager.h"
#include <cstring>

// Atom Parameter Manager (hash, atomParameterManager_enter, atomParameterManager_find)

AtomParameterManager::AtomParameterManager()
{
    memset(dictionary, 0, sizeof(dictionary));
}

AtomParameterManager::~AtomParameterManager()
{
    for (int i = 0; i < MAXKEY; i++)
        if (dictionary[i])
            delete dictionary[i];
}

unsigned int AtomParameterManager::hash(const char *key)
{
    switch (strlen(key))
    {
    case 0:
        return 0;
    case 1:
        return key[0];
    default:
        return (unsigned int)key[0] + 256*(unsigned int)key[1];
    }
}

void AtomParameterManager::insert(const char *key, const ParameterEntry &value)
{
    unsigned int hashKey = hash(key);
    if (!dictionary[hashKey])
        dictionary[hashKey] = new ParameterEntry();
    *dictionary[hashKey] = value;
}

ParameterEntry *AtomParameterManager::find(const char *key) const
{
    return dictionary[hash(key)];
}

int AtomParameterManager::getRecIndex(const char *key) const
{
    ParameterEntry *foundParam = find(key);
    if (foundParam)
        return foundParam->recIndex;
    return -1;
}
